# --- 0. Load Libraries & Setup ---
message("--- Loading Libraries ---")
suppressPackageStartupMessages({
    library(Seurat)
    library(mistyR)
    library(dplyr)
    library(tidyr)
    library(distances)
    library(ggplot2)
    library(future)
    library(patchwork)
    library(Matrix)
    library(SpatialExperiment)
    library(SingleCellExperiment)
    library(spacexr) # For RCTD
    library(tidyverse)
    # Optional libraries will be checked within functions if needed
})

# --- Helper Function: Save Plot ---
save_plot_pdf <- function(plot_object, filename, width = 7, height = 6, ...) {
    if (!is.null(plot_object) && inherits(plot_object, "ggplot")) {
        tryCatch({
            ggsave(filename = filename, plot = plot_object, width = width, height = height, device = "pdf", ...)
            message(paste("Saved plot:", filename))
            return(filename)
        }, error = function(e) {
            warning(paste("Could not save plot", filename, ":", e$message), immediate. = TRUE)
            return(NULL)
        })
    } else {
        warning(paste("Invalid plot object provided for:", filename), immediate. = TRUE)
        return(NULL)
    }
}


# --- MODULE 1: Seurat Preprocessing ---

#' Preprocess Spatial Seurat Object
#'
#' Loads 10X spatial data and performs standard Seurat preprocessing.
#' @param data_dir Path to 10X data directory.
#' @param rna_assay Name for the initial RNA assay.
#' @param sct_assay Name for the SCT assay.
#' @return Processed Seurat object.
preprocess_seurat_spatial <- function(data_dir, rna_assay = "Spatial", sct_assay = "SCT") {
    message("--- Loading Spatial Data & Preprocessing ---")
    h5_file <- file.path(data_dir, 'filtered_feature_bc_matrix.h5')
    if (!file.exists(h5_file)) stop(paste("H5 file not found at:", h5_file))

    tryCatch({
        seu_obj <- Load10X_Spatial(data.dir = data_dir, filename = basename(h5_file), assay = rna_assay)
    }, error = function(e) stop(paste("Error loading spatial data:", e$message)))

    tryCatch({
        n_rna_features <- 0
        if(rna_assay %in% names(seu_obj@assays)) n_rna_features <- nrow(seu_obj[[rna_assay]])

        if (sct_assay != rna_assay && n_rna_features > 0) {
             message("Running SCTransform...")
             seu_obj <- SCTransform(seu_obj, assay = rna_assay, verbose = FALSE)
        } else if (!sct_assay %in% names(seu_obj@assays) && n_rna_features > 0) {
             message("SCT assay not found, running NormalizeData/FindVariableFeatures/ScaleData on RNA assay...")
             seu_obj <- NormalizeData(seu_obj, assay = rna_assay, verbose = FALSE)
             seu_obj <- FindVariableFeatures(seu_obj, assay = rna_assay, verbose = FALSE)
             seu_obj <- ScaleData(seu_obj, assay = rna_assay, verbose = FALSE)
             sct_assay <- rna_assay # Use RNA assay name if SCT wasn't run
        } else if (!sct_assay %in% names(seu_obj@assays) && n_rna_features == 0) {
             stop("No counts found in the initial RNA assay.")
        }

        pca_features <- intersect(rownames(seu_obj), VariableFeatures(seu_obj[[sct_assay]]))
        if(length(pca_features) < 50) pca_features <- rownames(seu_obj[[sct_assay]])

        seu_obj <- RunPCA(seu_obj, assay = sct_assay, npcs = 30, features = pca_features, verbose = FALSE)
        seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
        seu_obj <- FindNeighbors(seu_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
        seu_obj <- FindClusters(seu_obj, resolution = 0.8, verbose = FALSE)
    }, error = function(e) stop(paste("Error during Seurat preprocessing:", e$message)))
    message("Seurat preprocessing complete.")
    return(seu_obj)
}


# --- MODULE 2: RCTD Deconvolution ---

#' Run RCTD Deconvolution
#'
#' Performs RCTD deconvolution using a single-cell reference.
#' @param seu_obj Processed spatial Seurat object.
#' @param ref_path Path to the reference Seurat object RDS file.
#' @param coords_df A data frame with spatial coordinates (barcode, row, col).
#' @param rna_assay Name of the assay holding raw counts in `seu_obj`.
#' @param ref_ident_col Column name for cell types in reference metadata.
#' @param ref_counts_col Column name for counts/UMI in reference metadata.
#' @param rctd_cores Number of cores for RCTD.
#' @return Seurat object with RCTD results added to metadata.
run_rctd_spatial <- function(seu_obj, ref_path, coords_df,
                             rna_assay = "Spatial",
                             ref_ident_col = "updated_celltype",
                             ref_counts_col = "nCount_RNA",
                             rctd_cores = 4) {
    message("--- Loading Reference & Running RCTD ---")
    if (!file.exists(ref_path)) stop(paste("Reference RDS file not found at:", ref_path))
    tryCatch({
        ref_seu <- readRDS(ref_path)
        if (!inherits(ref_seu, "Seurat")) {
             if (inherits(ref_seu, "SingleCellExperiment")) { ref_seu <- as.Seurat(ref_seu, counts = "counts", data = "logcounts") }
             else { stop("Reference file is not a Seurat object or SingleCellExperiment.") }
        }
        if (!"RNA" %in% names(ref_seu@assays)) stop("Reference must contain an 'RNA' assay.")
        if (!ref_ident_col %in% colnames(ref_seu@meta.data)) stop(paste("Reference metadata missing column:", ref_ident_col))
        if (!ref_counts_col %in% colnames(ref_seu@meta.data)) stop(paste("Reference metadata missing column:", ref_counts_col))
        ref_counts <- ref_seu[['RNA']]$counts
        ref_cluster <- as.factor(ref_seu@meta.data[[ref_ident_col]])
        names(ref_cluster) <- colnames(ref_seu)
        ref_nUMI <- ref_seu@meta.data[[ref_counts_col]]
        names(ref_nUMI) <- colnames(ref_seu)
        reference <- Reference(ref_counts, ref_cluster, ref_nUMI)
        rm(ref_seu, ref_counts, ref_cluster, ref_nUMI); gc()
    }, error = function(e) stop(paste("Error preparing RCTD reference:", e$message)))
    message("RCTD reference prepared.")

    tryCatch({
        query_counts <- seu_obj[[rna_assay]]$counts
        query_coords_prep <- coords_df %>% filter(barcode %in% colnames(query_counts)) %>% column_to_rownames("barcode")
        query_coords_prep <- query_coords_prep[colnames(query_counts), 1:2]
        colnames(query_coords_prep) <- c("x", "y")
        query <- SpatialRNA(query_coords_prep, query_counts, colSums(query_counts))
        RCTD <- create.RCTD(query, reference = reference, max_cores = rctd_cores, CELL_MIN_INSTANCE = 10)
        RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')
        rctd_results_df <- RCTD@results$results_df
        common_barcodes_rctd <- intersect(rownames(seu_obj@meta.data), rownames(rctd_results_df))
        saveRDS(list(RCTD,rctd_results_df),'RCTD_resultsList.rds')
        seu_obj <- AddMetaData(seu_obj, metadata = rctd_results_df[common_barcodes_rctd, , drop = FALSE])
        rm(query, reference, RCTD, rctd_results_df, query_coords_prep); gc()
    }, error = function(e) {
        warning(paste("Error during RCTD run:", e$message), immediate. = TRUE)
        if (!"first_type" %in% colnames(seu_obj@meta.data)) stop("RCTD failed or results not added.")
    })
    message("RCTD analysis complete and results added.")
    return(seu_obj)
}


# --- [Optional] MODULE 3: CellChat Analysis ---

#' Run CellChat Analysis on Spatial Data
#'
#' Performs CellChat analysis using spot expression and cell type labels.
#' @param seu_obj Seurat object.
#' @param coords_df A data frame with spatial coordinates (barcode, row, col).
#' @param rna_assay Assay containing log-normalized data for CellChat.
#' @param cell_identity_col Metadata column with cell type label per spot (e.g., 'first_type').
#' @param db CellChat database (e.g., CellChatDB.human).
#' @param keywords_pattern Regex pattern to filter pathways.
#' @return A tibble with CellChat pathway scores per barcode (barcode, CC_Pathway1, CC_Pathway2, ...). Returns empty tibble with barcode column on failure or if no pathways found.
run_cellchat_spatial <- function(seu_obj, coords_df,
                                 rna_assay = "RNA",
                                 cell_identity_col = "second_type",
                                 db = NULL,
                                 keywords_pattern = "TNF|IL[0-9]|IFN|TGF|CXCL|CCL|EDN|ANG|ADREN|GLUCAGON|CAMP|VASO") {
    message("--- [Optional] Starting CellChat Analysis ---")
    if (!requireNamespace("CellChat", quietly = TRUE)) {
        warning("CellChat package not found. Skipping CellChat analysis.", immediate. = TRUE)
        return(tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj))))
    }
    if (is.null(db)) {
        warning("CellChat database not provided (e.g., CellChatDB.human). Skipping CellChat analysis.", immediate. = TRUE)
        return(tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj))))
    }

    cellchat_scores_df <- tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj))) # Default empty

    tryCatch({
        expression_matrix_cc <- GetAssayData(seu_obj, assay = rna_assay, slot = "data")
        if (!cell_identity_col %in% colnames(seu_obj@meta.data)) {
             stop(paste("Cell identity column", cell_identity_col, "not found in Seurat metadata."))
        }
        cell_labels_cc <- seu_obj@meta.data[[cell_identity_col]]
        names(cell_labels_cc) <- rownames(seu_obj@meta.data)
        common_cells_cc <- intersect(colnames(expression_matrix_cc), names(cell_labels_cc))
        common_cells_cc <- intersect(coords_df$barcode, common_cells_cc)
        expression_matrix_cc_filt <- expression_matrix_cc[, common_cells_cc]
        cell_labels_cc_filt <- cell_labels_cc[common_cells_cc]
        stopifnot(ncol(expression_matrix_cc_filt) == length(cell_labels_cc_filt))
        stopifnot(all(colnames(expression_matrix_cc_filt) == names(cell_labels_cc_filt)))
        meta_df_cellchat <- data.frame(labels = cell_labels_cc_filt, row.names = names(cell_labels_cc_filt))
        meta_df_cellchat$labels = droplevels(meta_df_cellchat$labels, exclude = setdiff(levels(meta_df_cellchat$labels),unique(meta_df_cellchat$labels)))

        cellchat <- createCellChat(object = expression_matrix_cc_filt, meta = meta_df_cellchat, group.by = "labels")
        cellchat@DB <- db
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)
        saveRDS(cellchat,'spatial_cellchat_tmp.rds')

        pathways_all <- cellchat@netP$pathways
        pathways_interest <- grep(keywords_pattern, pathways_all, ignore.case = TRUE, value = TRUE)
        if (length(pathways_interest) > 0) {
            message(paste("Found relevant CellChat pathways:", paste(pathways_interest, collapse=", ")))
            prob_array <- cellchat@netP$prob
            prob_array_interest <- prob_array[,, pathways_interest, drop=FALSE]
            cell_group_levels <- levels(cellchat@idents)
            spot_scores_mat <- matrix(0, nrow = length(common_cells_cc), ncol = length(pathways_interest), dimnames = list(common_cells_cc, pathways_interest))
            target_type_indices <- match(cell_labels_cc_filt, cell_group_levels)
            valid_indices <- !is.na(target_type_indices)
            target_type_indices_valid <- target_type_indices[valid_indices]

            for (pathway_idx in 1:length(pathways_interest)) {
                prob_sum_per_target_type <- apply(prob_array_interest[,, pathway_idx, drop = FALSE], MARGIN = 2, FUN = sum)
                # Ensure prob_sum_per_target_type has names matching cell_group_levels if possible
                if(!is.null(names(prob_sum_per_target_type)) && all(cell_group_levels %in% names(prob_sum_per_target_type))){
                   spot_scores_mat[valid_indices, pathway_idx] <- prob_sum_per_target_type[cell_group_levels[target_type_indices_valid]]
                } else {
                   # Fallback if names are missing or don't match - relies on order
                   warning("Target type names missing or mismatched in CellChat probability sum. Assigning scores based on index order.", immediate. = TRUE)
                   spot_scores_mat[valid_indices, pathway_idx] <- prob_sum_per_target_type[target_type_indices_valid]
                }
            }
            cellchat_scores_df <- as_tibble(spot_scores_mat, rownames = "barcode") %>%
              rename_with(~paste0("CC_", str_replace_all(., "[^[:alnum:]_]", "_")), -barcode)
            message("Summary of extracted CellChat scores (first few pathways):")
            print(summary(cellchat_scores_df[, 2:min(6, ncol(cellchat_scores_df)), drop=FALSE]))
        } else {
            warning("No pathways matching the keywords found in CellChat results.", immediate. = TRUE)
        }
        rm(cellchat, expression_matrix_cc, cell_labels_cc, expression_matrix_cc_filt, cell_labels_cc_filt, meta_df_cellchat, prob_array, prob_array_interest, spot_scores_mat); gc()
    }, error = function(e) {
        warning(paste("Error during CellChat analysis:", e$message), immediate. = TRUE)
    })
    message("--- Finished CellChat Analysis ---")
    # Ensure barcode column exists even if analysis failed partially
    if (!"barcode" %in% colnames(cellchat_scores_df)) {
        cellchat_scores_df <- tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj)))
    }
    return(cellchat_scores_df)
}


# --- [Optional] MODULE 4: Dorothea/Viper Analysis ---

#' Run Dorothea/Viper Analysis on Spatial Data
#'
#' Calculates TF activity using Viper and Dorothea regulons.
#' @param seu_obj Seurat object.
#' @param coords_df A data frame with spatial coordinates (barcode, row, col).
#' @param rna_assay Assay containing appropriate expression data for Viper (e.g., log-normalized).
#' @param regulon Dorothea regulon object (e.g., filtered `dorothea_hs`).
#' @param tfs_list A list containing named vectors of TFs to extract.
#' @return A tibble with TF activity scores per barcode (barcode, TF_TF1, TF_TF2, ...). Returns empty tibble with barcode column on failure or if no TFs found.
run_viper_spatial <- function(seu_obj, coords_df,
                              rna_assay = "RNA",
                              regulon = NULL,
                              tfs_list = list(
                                  inflammation=c("NFKB1", "RELA", "STAT1", "STAT3", "IRF1", "IRF7", "JUN", "FOS"),
                                  vasoconstriction=c("GATA2", "GATA3", "SRF", "MEF2C", "ETV1"),
                                  cAMP=c("CREB1", "ATF1", "ATF4")
                              )) {
    message("--- [Optional] Starting Dorothea/Viper Analysis ---")
    library('dorothea')
    if (!requireNamespace("dorothea", quietly = TRUE) || !requireNamespace("viper", quietly = TRUE)) {
        warning("Dorothea or Viper package not found. Skipping Dorothea/Viper analysis.", immediate. = TRUE)
        return(tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj))))
    }
     if (is.null(regulon)) {
        warning("Dorothea regulon not provided. Loading default human A/B/C.", immediate. = TRUE)
        dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
        regulon <- dorothea_regulon_human %>% filter(confidence %in% c("A", "B", "C"))
     }
   

    dorothea_scores_df <- tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj))) # Default empty

    tryCatch({
        expression_matrix_viper <- GetAssayData(seu_obj, assay = rna_assay, slot = "data")
        common_cells_viper <- intersect(colnames(expression_matrix_viper), coords_df$barcode)
        expression_matrix_viper_filt <- expression_matrix_viper[, common_cells_viper]

        tf_activity <- run_viper(as.matrix(expression_matrix_viper_filt), regulon, options =  list(method = "scale", minsize = 4,
                                                                                        eset.filter = FALSE, cores = 1,
                                                                                        verbose = FALSE))

        tfs_interest_names <- unique(unlist(tfs_list))
        tfs_available <- intersect(tfs_interest_names, rownames(tf_activity))
        if (length(tfs_available) > 0) {
            message(paste("Found relevant TFs:", paste(tfs_available, collapse=", ")))
            dorothea_scores <- t(tf_activity[tfs_available, , drop = FALSE])
            dorothea_scores_df <- dorothea_scores %>%
              as_tibble(rownames = "barcode") %>%
              rename_with(~paste0("TF_", .), -barcode) %>%
              filter(barcode %in% common_cells_viper)
            message("Summary of extracted Dorothea TF activity scores (first few TFs):")
            print(summary(dorothea_scores_df[, 2:min(6, ncol(dorothea_scores_df)), drop=FALSE]))
        } else {
            warning("None of the specified TFs of interest were found in Dorothea/Viper results.", immediate. = TRUE)
        }
        rm(expression_matrix_viper, expression_matrix_viper_filt, tf_activity, dorothea_scores); gc()
    }, error = function(e) {
        warning(paste("Error during Dorothea/Viper analysis:", e$message), immediate. = TRUE)
    })
    message("--- Finished Dorothea/Viper Analysis ---")
    # Ensure barcode column exists even if analysis failed partially
    if (!"barcode" %in% colnames(dorothea_scores_df)) {
        dorothea_scores_df <- tibble(barcode = intersect(coords_df$barcode, colnames(seu_obj)))
    }
    return(dorothea_scores_df)
}


# --- [Optional] MODULE 5: Banksy Analysis ---

#' Run Banksy Analysis on Spatial Data
#'
#' Performs Banksy domain identification using HVGs.
#' @param seu_obj Seurat object (SCT assay recommended).
#' @param coords_df A data frame with spatial coordinates (barcode, row, col).
#' @param sct_assay Assay containing normalized data (e.g., SCT data slot).
#' @param k_geom K for Banksy spatial graph.
#' @param lambda Lambda for Banksy.
#' @param k_centers K for Banksy clustering.
#' @param output_dir Directory to save Banksy plot.
#' @param sample_name Sample name for plot filename.
#' @return Seurat object with 'banksy_domain' metadata added. Returns original object if Banksy fails.
run_banksy_spatial <- function(seu_obj, coords_df,
                               sct_assay = "SCT",
                               k_geom = 15, lambda = 0.2, k_centers = 10,
                               output_dir = ".", sample_name = "sample") {
    message("--- [Optional] Starting Banksy Analysis ---")
     if (!requireNamespace("Banksy", quietly = TRUE) || !requireNamespace("SpatialExperiment", quietly = TRUE)) {
        warning("Banksy or SpatialExperiment package not found. Skipping Banksy analysis.", immediate. = TRUE)
        return(seu_obj)
    }

    original_seu_obj <- seu_obj # Keep original in case of error
    tryCatch({
        message("Finding variable features for Banksy...")
        if (!sct_assay %in% names(seu_obj@assays)) stop(paste("SCT assay", sct_assay, "not found for Banksy HVG calculation."))
        if (length(VariableFeatures(seu_obj[[sct_assay]])) < 100) {
             seu_obj <- FindVariableFeatures(seu_obj, assay = sct_assay, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
        }
        hvgs <- VariableFeatures(seu_obj, assay = sct_assay)
        hvg_expr_matrix <- GetAssayData(seu_obj, assay = sct_assay, layer = "data")[hvgs, ] %>% as.matrix()
        common_barcodes_banksy <- intersect(colnames(hvg_expr_matrix), coords_df$barcode)
        hvg_expr_matrix <- hvg_expr_matrix[, common_barcodes_banksy]
        banksy_coords <- coords_df %>% filter(barcode %in% common_barcodes_banksy)
        stopifnot(ncol(hvg_expr_matrix) == nrow(banksy_coords), all(colnames(hvg_expr_matrix) == banksy_coords$barcode))

        spatial_coords_banksy <- as.matrix(banksy_coords %>% select(row, col))
        rownames(spatial_coords_banksy) <- banksy_coords$barcode
        spe <- SpatialExperiment(assays = list(logcounts = hvg_expr_matrix), spatialCoords = spatial_coords_banksy)

        message("Running computeBanksy...")
        spe <- computeBanksy(spe, assay_name = "logcounts", M = 1, k_geom = k_geom, lambda = lambda, verbose = FALSE)
        message("Running runBanksyPCA...")
        spe <- runBanksyPCA(spe, M = 1, lambda = lambda, npcs = 20)
        message("Running runBanksyUMAP...")
        spe <- runBanksyUMAP(spe, M = 1, lambda = lambda, npcs = 20)
        message("Running clusterBanksy...")
        spe <- clusterBanksy(spe, M = 1, lambda = lambda, resolution = 1.0, k_centers = k_centers)

        cluster_col_name <- grep(paste0("clust_M1_lam", lambda), colnames(colData(spe)), value = TRUE)
        if (length(cluster_col_name) == 0) cluster_col_name <- grep(paste0("clust_M1_lam", lambda, "_k", k_centers), colnames(colData(spe)), value = TRUE)
        if (length(cluster_col_name) == 0) stop("Could not find Banksy cluster results column.")
        cluster_col_name <- cluster_col_name[1]
        message(paste("Using Banksy cluster column:", cluster_col_name))
        banksy_clusters <- colData(spe)[[cluster_col_name]]
        names(banksy_clusters) <- colnames(spe)

        # Visualize
        banksy_plot_data <- banksy_coords %>%
          inner_join(tibble(barcode = names(banksy_clusters), banksy_domain = factor(banksy_clusters)), by = "barcode")
        p_banksy_domains <- ggplot(banksy_plot_data, aes(x = col, y = row, color = banksy_domain)) +
          geom_point(size = 1.5, alpha = 0.8) + scale_y_reverse() + coord_fixed() +
          labs(title = paste(sample_name, "- Banksy Domains (HVG Expr)"), color = "Domain") + theme_minimal() +
          theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "bottom")
        banksy_plot_path <- file.path(output_dir, paste0(sample_name, "_banksy_domains.pdf")) # Save as PDF
        save_plot_pdf(p_banksy_domains, banksy_plot_path, width = 7, height = 6)

        # Add to Seurat object
        common_barcodes_seu <- intersect(rownames(seu_obj@meta.data), names(banksy_clusters))
        seu_obj <- AddMetaData(seu_obj, metadata = banksy_clusters[common_barcodes_seu], col.name = "banksy_domain")
        message("Banksy domains added to seu_obj@meta.data$banksy_domain")
        rm(spe, hvgs, hvg_expr_matrix, banksy_coords, spatial_coords_banksy, banksy_clusters, banksy_plot_data, p_banksy_domains); gc()
    }, error = function(e) {
        warning(paste("Error during Banksy analysis:", e$message), immediate. = TRUE)
        return(original_seu_obj) # Return original object on error
    })
    message("--- Finished Banksy Analysis ---")
    return(seu_obj)
}


# --- MODULE 6: Prepare mistyR Input Data ---

#' Prepare Input Data for mistyR Analysis
#'
#' Combines coordinates, target gene expression, cell type indicators,
#' and optional scores, handles NAs, and filters low variance features.
#'
#' @param seu_obj Seurat object (used for filtering barcodes if needed).
#' @param coords_df Data frame with coordinates (barcode, row, col).
#' @param target_gene The target gene name.
#' @param sct_assay Assay containing normalized target gene data.
#' @param cell_type_indicators Tibble with barcode and is_Myeloid, is_Arterial, is_Venous columns.
#' @param cellchat_scores_df Tibble with barcode and CellChat scores (starting with CC_). Can be NULL or empty.
#' @param dorothea_scores_df Tibble with barcode and Dorothea scores (starting with TF_). Can be NULL or empty.
#' @param min_spots Minimum number of spots required after NA removal.
#' @return A list containing `mistyR_pos` (positions) and `mistyR_expr` (features), or NULL if insufficient spots.
prepare_mistyR_input <- function(seu_obj, coords_df, target_gene, sct_assay,
                                 cell_type_indicators,
                                 cellchat_scores_df = NULL,
                                 dorothea_scores_df = NULL,
                                 min_spots = 50) {

    message("Preparing mistyR input data...")
    # Get target gene expression
    if (!target_gene %in% rownames(seu_obj[[sct_assay]])) {
        stop(paste("Target gene '", target_gene, "' not found in the", sct_assay, "assay."))
    }
    target_gene_expr_misty <- FetchData(seu_obj, vars = target_gene, assay = sct_assay, layer = 'data') %>%
                              as_tibble(rownames = "barcode")

    # Start building the input data frame
    mistyR_input_data <- inner_join(coords_df, target_gene_expr_misty, by = "barcode") %>%
                         inner_join(cell_type_indicators, by = "barcode")

    # Conditionally join CellChat scores
    if (!is.null(cellchat_scores_df) && nrow(cellchat_scores_df) > 0 && ncol(cellchat_scores_df) > 1) {
        mistyR_input_data <- left_join(mistyR_input_data, cellchat_scores_df, by = "barcode")
    }
    # Conditionally join Dorothea scores
    if (!is.null(dorothea_scores_df) && nrow(dorothea_scores_df) > 0 && ncol(dorothea_scores_df) > 1) {
        mistyR_input_data <- left_join(mistyR_input_data, dorothea_scores_df, by = "barcode")
    }

    # Handle NAs and check spot count
    mistyR_input_data <- mistyR_input_data %>% drop_na()
    if (nrow(mistyR_input_data) < min_spots) {
        warning(paste("Not enough spots remaining after merging/NA removal for mistyR:", nrow(mistyR_input_data)), immediate. = TRUE)
        return(NULL)
    }
    message(paste("Number of spots for mistyR:", nrow(mistyR_input_data)))

    # Separate positions and features
    mistyR_pos <- mistyR_input_data %>% select(row, col)
    mistyR_expr <- mistyR_input_data %>%
        select(
            any_of(target_gene),
            starts_with("is_"),
            starts_with("CC_"),
            starts_with("TF_")
        ) %>%
        mutate(across(everything(), as.numeric))

    # Check variance
    if(ncol(mistyR_expr) > 1) {
        variances <- apply(mistyR_expr, 2, var, na.rm = TRUE)
        zero_var_cols <- names(variances[variances == 0 | is.na(variances)])
        if(length(zero_var_cols) > 0) {
            warning(paste("Removing zero/NA variance columns:", paste(zero_var_cols, collapse=", ")), immediate. = TRUE)
            mistyR_expr <- mistyR_expr %>% select(-all_of(zero_var_cols))
        }
    }
    if(ncol(mistyR_expr) <= 1) {
        warning(paste("Not enough features remaining after variance check:", ncol(mistyR_expr)), immediate. = TRUE)
        return(NULL)
    }
    message(paste("Number of features for mistyR:", ncol(mistyR_expr)))
    # message("Final features for mistyR:"); print(colnames(mistyR_expr))

    return(list(mistyR_pos = mistyR_pos, mistyR_expr = mistyR_expr))
}


# --- MODULE 7: Run mistyR Core Analysis ---

#' Run Core mistyR Workflow
#'
#' Creates views, runs misty model training, and collects results.
#' @param mistyR_pos Data frame of spatial positions (row, col).
#' @param mistyR_expr Data frame of features for mistyR.
#' @param l_param Neighborhood radius for paraview.
#' @param output_folder Path to save raw mistyR results.
#' @param seed Random seed for reproducibility.
#' @return misty.results object, or NULL on failure.
run_mistyR_core <- function(mistyR_pos, mistyR_expr, l_param = 10, output_folder = "mistyR_run", seed = 123) {
    misty_results_obj <- NULL
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    tryCatch({
        misty.intra <- create_initial_view(mistyR_expr)
        misty.views <- misty.intra %>% add_paraview(mistyR_pos, l = l_param)
        message("mistyR views created.")
        message(paste("Running mistyR, results will be saved in:", output_folder))
        run_misty(misty.views, results.folder = output_folder, seed = seed)
        message("Collecting mistyR results...")
        misty_results_obj <- collect_results(output_folder)
        message("mistyR results collected.")
    }, error = function(e) {
        warning(paste("Error during mistyR analysis run:", e$message), immediate. = TRUE)
        misty_results_obj <<- NULL # Ensure return is NULL on error
    })
    return(misty_results_obj)
}


# --- MODULE 8: Plot mistyR Results ---

#' Generate and Save Standard mistyR Plots
#'
#' Creates contribution and heatmap plots from mistyR results.
#' @param misty_results The results object from `collect_results`.
#' @param plot_subdir Directory to save the plots.
#' @param plot_prefix Prefix for plot filenames (e.g., sample name or sample_domain).
#' @param l_param Neighborhood radius used (for titles/filenames).
#' @param cutoff_intra Cutoff for intraview heatmap.
#' @param cutoff_para Cutoff for paraview/contrast heatmaps.
#' @param file_format "pdf" or "png".
#' @return A list of paths to the saved plot files.
plot_mistyR_results <- function(misty_results, plot_subdir = ".", plot_prefix = "mistyR",
                                l_param = 10, cutoff_intra = 0.5, cutoff_para = 0.3,
                                file_format = "pdf") {

    if (is.null(misty_results) || length(misty_results) == 0) {
        warning("No mistyR results provided for plotting.", immediate. = TRUE)
        return(list())
    }
    if (!dir.exists(plot_subdir)) dir.create(plot_subdir, recursive = TRUE)

    message(paste("Generating mistyR plots for prefix:", plot_prefix))
    plots_saved <- list()
    file_ext <- paste0(".", tolower(file_format))

    # Define plot paths
    plot_path_contrib <- file.path(plot_subdir, paste0(plot_prefix, "_contribution", file_ext))
    plot_path_intra <- file.path(plot_subdir, paste0(plot_prefix, "_heatmap_intra", file_ext))
    plot_path_para <- file.path(plot_subdir, paste0(plot_prefix, "_heatmap_para", file_ext))
    plot_path_contrast <- file.path(plot_subdir, paste0(plot_prefix, "_heatmap_contrast", file_ext))

    # Estimate plot parameters based on number of features
    n_features_plot <- tryCatch(ncol(misty_results$importances.aggregated$importance), error = function(e) 4) # Default if error
    heatmap_dims <- max(6, n_features_plot * 0.4)
    heatmap_text_size <- max(3, 8 - floor(n_features_plot/10))

    # Contribution Plot
    tryCatch({
        p_contrib <- misty_results %>% plot_view_contributions() + labs(title = paste(plot_prefix, "- View Contributions")) + theme(plot.title = element_text(size=10))
        saved_path <- save_plot_pdf(p_contrib, plot_path_contrib, width = 6, height = 5, device = ifelse(file_format=="pdf", "pdf", "png"), bg="white")
        if(!is.null(saved_path)) plots_saved[["contribution"]] <- saved_path
    }, error = function(e) { warning(paste("Could not generate/save contribution plot:", e$message), immediate. = TRUE) })

    # Heatmaps
    para_view_name_plot <- names(misty_results)[grep("^para", names(misty_results))][1]
    if (is.na(para_view_name_plot)) {
        warning("Paraview results not found in misty_results object for heatmaps.", immediate. = TRUE)
    } else {
        # Intraview
        tryCatch({
            p_heatmap_intra <- misty_results %>% plot_interaction_heatmap(view = "intra", cutoff = cutoff_intra)
            if (inherits(p_heatmap_intra, "ggplot")) {
                p_heatmap_intra <- p_heatmap_intra + labs(title = paste(plot_prefix, "- Intraview (Cutoff=", cutoff_intra, ")")) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, size=heatmap_text_size), axis.text.y = element_text(size=heatmap_text_size))
                saved_path <- save_plot_pdf(p_heatmap_intra, plot_path_intra, width = heatmap_dims, height = heatmap_dims*0.8, device = ifelse(file_format=="pdf", "pdf", "png"), bg="white")
                if(!is.null(saved_path)) plots_saved[["intraview"]] <- saved_path
            } else { warning("Intraview plot could not be generated.", immediate. = TRUE) }
        }, error = function(e) { warning(paste("Error generating/saving intraview heatmap:", e$message), immediate. = TRUE) })
        # Paraview
        tryCatch({
            p_heatmap_para <- misty_results %>% plot_interaction_heatmap(view = "para.10", cutoff = cutoff_para)
            if (inherits(p_heatmap_para, "ggplot")) {
                p_heatmap_para <- p_heatmap_para + labs(title = paste(plot_prefix, "- Paraview (l=", l_param, ", Cutoff=", cutoff_para, ")")) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, size=heatmap_text_size), axis.text.y = element_text(size=heatmap_text_size))
                saved_path <- save_plot_pdf(p_heatmap_para, plot_path_para, width = heatmap_dims, height = heatmap_dims*0.8, device = ifelse(file_format=="pdf", "pdf", "png"), bg="white")
                if(!is.null(saved_path)) plots_saved[["paraview"]] <- saved_path
            } else { warning("Paraview plot could not be generated.", immediate. = TRUE) }
        }, error = function(e) { warning(paste("Error generating/saving paraview heatmap:", e$message), immediate. = TRUE) })
        # Contrast
        tryCatch({
            p_heatmap_contrast <- misty_results %>% plot_contrast_heatmap("intra", "para.10", cutoff = cutoff_para)
            if (inherits(p_heatmap_contrast, "ggplot")) {
                 p_heatmap_contrast <- p_heatmap_contrast + labs(title = paste(plot_prefix, "- Contrast (Cutoff=", cutoff_para, ")")) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, size=heatmap_text_size), axis.text.y = element_text(size=heatmap_text_size))
                saved_path <- save_plot_pdf(p_heatmap_contrast, plot_path_contrast, width = heatmap_dims, height = heatmap_dims*0.8, device = ifelse(file_format=="pdf", "pdf", "png"), bg="white")
                if(!is.null(saved_path)) plots_saved[["contrast"]] <- saved_path
            } else { warning("Contrast plot could not be generated.", immediate. = TRUE) }
        }, error = function(e) { warning(paste("Error generating/saving contrast heatmap:", e$message), immediate. = TRUE) })
    }
    return(plots_saved)
}


# --- Main Orchestrating Function ---

#' Run Full Spatial Analysis Workflow (Modular)
#'
#' Main function to orchestrate the spatial analysis workflow, calling modular
#' functions for preprocessing, deconvolution, optional analyses (Banksy, CellChat,
#' Dorothea), and mistyR.
#'
#' @inheritParams preprocess_seurat_spatial
#' @inheritParams run_rctd_spatial
#' @inheritParams run_cellchat_spatial
#' @inheritParams run_viper_spatial
#' @inheritParams run_banksy_spatial
#' @inheritParams prepare_mistyR_input
#' @inheritParams run_mistyR_core
#' @inheritParams plot_mistyR_results
#' @param output_dir Path for all outputs.
#' @param sample_name Sample identifier.
#' @param target_gene Target gene for mistyR.
#' @param cell_type_mapping List mapping generic to specific cell type names.
#' @param run_banksy Logical, run Banksy?
#' @param run_cellchat Logical, run CellChat?
#' @param run_dorothea Logical, run Dorothea/Viper?
#' @param mistyR_per_domain Logical, run mistyR per Banksy domain?
#' @param n_cores Cores for future parallelization.
#' @param plot_file_format "pdf" or "png" for mistyR plots.
#'
#' @return Path to the final saved Seurat object RDS file.
run_spatial_analysis_modular <- function(
    data_dir, ref_path, output_dir, sample_name,
    target_gene = "RAPGEF3",
    cell_type_mapping = list(
        Myeloid = "Immune_Myeloid", Arterial_Pos = "Arterial_Endo_RAPGEF3+SELP+",
        Arterial_Neg = "Arterial_RAPGEF3-_Cells", Venous_Pos = "Venous_Endo_RAPGEF3+SELP-",
        Venous_Neg = "Venous_RAPGEF3-_Cells"
    ),
    run_banksy = FALSE, run_cellchat = FALSE, run_dorothea = FALSE, mistyR_per_domain = FALSE,
    cellchat_keywords_pattern = "TNF|IL[0-9]|IFN|TGF|CXCL|CCL|EDN|ANG|ADREN|GLUCAGON|CAMP|VASO",
    dorothea_tfs_list = list(
        inflammation=c("NFKB1", "RELA", "STAT1", "STAT3", "IRF1", "IRF7", "JUN", "FOS"),
        vasoconstriction=c("GATA2", "GATA3", "SRF", "MEF2C", "ETV1"),
        cAMP=c("CREB1", "ATF1", "ATF4")
    ),
    misty_l_param = 10, misty_cutoff_intra = 0.5, misty_cutoff_para = 0.3,
    n_cores = max(1, future::availableCores() - 2), rctd_cores = 4,
    sct_assay = "SCT", rna_assay = "Spatial",
    ref_ident_col = "updated_celltype", ref_counts_col = "nCount_RNA",
    banksy_k_geom = 15, banksy_lambda = 0.2, banksy_k_centers = 10,
    plot_file_format = "pdf"
    ) {

    # --- Setup ---
    message(paste("--- Starting Full Workflow for Sample:", sample_name, "---"))
    # Setup parallel backend
    plan(multisession, workers = n_cores)
    options(future.globals.maxSize = 8 * 1024^3)
    options(stringsAsFactors = FALSE)
    # Create output directory
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    # Validate flags
    if (mistyR_per_domain && !run_banksy) {
        warning("mistyR_per_domain=TRUE requires run_banksy=TRUE. Setting mistyR_per_domain=FALSE.", immediate. = TRUE)
        mistyR_per_domain <- FALSE
    }

    # --- 1. Preprocessing ---
    seu_obj <- preprocess_seurat_spatial(data_dir = data_dir, rna_assay = rna_assay, sct_assay = sct_assay)

    # --- 1.1 Coordinates ---
    coords_df <- NULL
    tryCatch({
        stopifnot(length(names(seu_obj@images)) > 0)
        image_key <- names(seu_obj@images)[1]
        coords_raw <- GetTissueCoordinates(seu_obj, image = image_key)
        if (is.null(coords_raw) || nrow(coords_raw) == 0) stop("GetTissueCoordinates failed.")
        coords_df <- coords_raw %>% as_tibble(rownames = "barcode")
        if (all(c("x", "y") %in% colnames(coords_df))) {
            coords_df <- coords_df %>% select(barcode, x, y)
            colnames(coords_df) <- c("barcode", "row", "col")
        } else if (all(c("imagerow", "imagecol") %in% colnames(coords_df))) {
            coords_df <- coords_df %>% select(barcode, imagerow, imagecol)
            colnames(coords_df) <- c("barcode", "row", "col")
        } else if (all(c("row", "col") %in% colnames(coords_df))) {
             if (!"barcode" %in% colnames(coords_df)) stop("Barcode missing.")
             coords_df <- coords_df %>% select(barcode, row, col)
        } else { stop("Standard coordinate columns not found.") }
        stopifnot(is.data.frame(coords_df), nrow(coords_df) > 0)
        message("'coords_df' created successfully.")
    }, error = function(e) stop(paste("Error getting coordinates:", e$message)))


    # --- 2. RCTD ---
    seu_obj <- run_rctd_spatial(seu_obj = seu_obj, ref_path = ref_path, coords_df = coords_df,
                                rna_assay = rna_assay, ref_ident_col = ref_ident_col,
                                ref_counts_col = ref_counts_col, rctd_cores = rctd_cores)

    # --- 3. CellChat (Optional) ---
    cellchat_scores_df <- NULL
    if(run_cellchat) {
        # Assuming human DB, adjust if needed
        library(CellChat)
        cellchat_scores_df <- run_cellchat_spatial(seu_obj = seu_obj, coords_df = coords_df,
                                                   rna_assay = sct_assay, # Use RNA assay for CellChat
                                                   cell_identity_col = 'second_type',
                                                   db = CellChatDB.human,
                                                   keywords_pattern = cellchat_keywords_pattern)
    }

    # --- 4. Dorothea/Viper (Optional) ---
    dorothea_scores_df <- NULL
    if(run_dorothea) {
        # Load default human regulon if not provided (adjust if needed)
         dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
         regulon_viper <- dorothea_regulon_human %>% filter(confidence %in% c("A", "B", "C"))
         dorothea_scores_df <- run_viper_spatial(seu_obj = seu_obj, coords_df = coords_df,
                                                rna_assay = rna_assay, # Use RNA assay for Viper
                                                regulon = regulon_viper,
                                                tfs_list = dorothea_tfs_list)
    }

    # --- 5. Banksy (Optional) ---
    if(run_banksy) {
        seu_obj <- run_banksy_spatial(seu_obj = seu_obj, coords_df = coords_df,
                                      sct_assay = sct_assay, # Use SCT assay for HVGs
                                      k_geom = banksy_k_geom, lambda = banksy_lambda, k_centers = banksy_k_centers,
                                      output_dir = output_dir, sample_name = sample_name)
    }

    # --- 6. Prepare mistyR Input ---
    # Create cell type indicators based on mapping
    rctd_types_misty <- seu_obj@meta.data %>% as_tibble(rownames = "barcode") %>% select(barcode, first_type, second_type)
    cell_type_indicators_misty <- rctd_types_misty %>%
      mutate(
        is_Myeloid = if_else(first_type == cell_type_mapping$Myeloid | second_type == cell_type_mapping$Myeloid, 1, 0),
        is_Arterial = if_else(first_type == cell_type_mapping$Arterial_Pos | second_type == cell_type_mapping$Arterial_Pos | first_type == cell_type_mapping$Arterial_Neg | second_type == cell_type_mapping$Arterial_Neg, 1, 0),
        is_Venous = if_else(first_type == cell_type_mapping$Venous_Pos | second_type == cell_type_mapping$Venous_Pos | first_type == cell_type_mapping$Venous_Neg | second_type == cell_type_mapping$Venous_Neg, 1, 0)
      ) %>% select(barcode, is_Myeloid, is_Arterial, is_Venous) # Adjust names based on mapping keys if needed

    mistyR_prepared_input <- prepare_mistyR_input(
        seu_obj = seu_obj, coords_df = coords_df, target_gene = target_gene, sct_assay = sct_assay,
        cell_type_indicators = cell_type_indicators_misty,
        cellchat_scores_df = if(run_cellchat) cellchat_scores_df else NULL,
        dorothea_scores_df = if(run_dorothea) dorothea_scores_df else NULL
    )

    # --- 7. Run mistyR (Global or Per Domain) ---
    mistyR_results <- NULL
    mistyR_plots_saved <- list()
    plot_subdir <- file.path(output_dir, paste0(sample_name, "_mistyR_plots"))

    if (!is.null(mistyR_prepared_input)) {
        if (mistyR_per_domain) {
            message("--- Running mistyR Per Domain ---")
            if (!"banksy_domain" %in% colnames(seu_obj@meta.data)) {
                stop("Cannot run mistyR per domain because Banksy analysis was skipped or failed.")
            }
            unique_domains <- na.omit(unique(seu_obj$banksy_domain))
            all_misty_results_list <- list()
            min_spots_misty_domain <- 50

            for (domain_label in unique_domains) {
                current_domain_label_str <- as.character(domain_label)
                message(paste("--- Processing Domain:", current_domain_label_str, "---"))
                barcodes_in_domain <- rownames(seu_obj@meta.data)[which(seu_obj$banksy_domain == domain_label)]
                if (length(barcodes_in_domain) < min_spots_misty_domain) {
                    message(paste("Skipping domain", current_domain_label_str, "- too few spots."))
                    next
                }
                # Prepare subset input data
                mistyR_input_subset <- prepare_mistyR_input(
                    seu_obj = seu_obj[, barcodes_in_domain], # Subset Seurat object for FetchData
                    coords_df = coords_df %>% filter(barcode %in% barcodes_in_domain),
                    target_gene = target_gene, sct_assay = sct_assay,
                    cell_type_indicators = cell_type_indicators_misty %>% filter(barcode %in% barcodes_in_domain),
                    cellchat_scores_df = if(run_cellchat) cellchat_scores_df %>% filter(barcode %in% barcodes_in_domain) else NULL,
                    dorothea_scores_df = if(run_dorothea) dorothea_scores_df %>% filter(barcode %in% barcodes_in_domain) else NULL,
                    min_spots = min_spots_misty_domain * 0.6
                 )

                 if (!is.null(mistyR_input_subset)) {
                    mistyR_domain_folder <- file.path(output_dir, paste0(sample_name, "_mistyR_per_domain"), paste0("domain_", current_domain_label_str))
                    domain_results <- run_mistyR_core(
                        mistyR_pos = mistyR_input_subset$mistyR_pos,
                        mistyR_expr = mistyR_input_subset$mistyR_expr,
                        l_param = misty_l_param,
                        output_folder = mistyR_domain_folder,
                        seed = 123
                    )
                    if (!is.null(domain_results)) {
                        all_misty_results_list[[paste0("domain_", current_domain_label_str)]] <- domain_results
                        # Plot results for this domain
                        domain_plots <- plot_mistyR_results(
                            misty_results = domain_results,
                            plot_subdir = plot_subdir,
                            plot_prefix = paste0(sample_name, "_domain_", current_domain_label_str),
                            l_param = misty_l_param, cutoff_intra = misty_cutoff_intra, cutoff_para = misty_cutoff_para,
                            file_format = plot_file_format
                         )
                         mistyR_plots_saved <- c(mistyR_plots_saved, domain_plots)
                    }
                 } else {
                    message(paste("Skipping domain", current_domain_label_str, "due to insufficient data after preparation."))
                 }
            } # End domain loop
            mistyR_results <- all_misty_results_list # Assign list for potential return/further use

        } else {
            # --- Run Global mistyR ---
            message("--- Running Global mistyR Analysis ---")
            mistyR_global_folder <- file.path(output_dir, paste0(sample_name, "_mistyR_global_results"))
            mistyR_results <- run_mistyR_core(
                 mistyR_pos = mistyR_prepared_input$mistyR_pos,
                 mistyR_expr = mistyR_prepared_input$mistyR_expr,
                 l_param = misty_l_param,
                 output_folder = mistyR_global_folder,
                 seed = 123
             )
             # Plot global results
             global_plots <- plot_mistyR_results(
                 misty_results = mistyR_results,
                 plot_subdir = plot_subdir,
                 plot_prefix = paste0(sample_name, "_global"),
                 l_param = misty_l_param, cutoff_intra = misty_cutoff_intra, cutoff_para = misty_cutoff_para,
                 file_format = plot_file_format
             )
             mistyR_plots_saved <- c(mistyR_plots_saved, global_plots)
        }
    } else {
         message("Skipping mistyR analysis because input data preparation failed.")
    }

    # --- 8. Final Save & Cleanup ---
    seurat_output_path <- file.path(output_dir, paste0(sample_name, "_processed_seu_final.rds"))
    message(paste("--- Saving Final Seurat Object to:", seurat_output_path, "---"))
    tryCatch({
        # Ensure optional metadata exists before saving if it was supposed to run
        if(run_banksy && !"banksy_domain" %in% colnames(seu_obj@meta.data)) warning("Banksy was set to run but 'banksy_domain' metadata is missing.", immediate.=TRUE)
        saveRDS(seu_obj, seurat_output_path)
    }, error = function(e) {
        warning(paste("Error saving final Seurat object:", e$message), immediate. = TRUE)
    })

    message("--- Cleaning up ---")
    plan(sequential)
    gc()

    message(paste("--- Workflow Complete for Sample:", sample_name, "---"))
    if(file.exists(seurat_output_path)) message(paste("Final Seurat object saved to:", seurat_output_path))
    if(length(mistyR_plots_saved) > 0) {
        message("mistyR plots saved to directory:", plot_subdir)
    } else {
        message("No mistyR plots were generated or saved.")
    }

    invisible(seurat_output_path)
}


# --- Example Usage ---
# # Define your parameters
# spatial_data_folder <- "FX0028/UA_HUTER_sp_10879894/"
# reference_rds <- "../sce_sub_for_spatial.rds"
# results_directory <- "Analysis_Results/FX0028_Modular"
# sample_id <- "FX0028_Modular"
# gene_of_interest <- "RAPGEF3"
# # IMPORTANT: Verify these names match your reference data's 'updated_celltype' column
cell_mapping <- list(
     Myeloid = "Immune_Myeloid",
     Arterial_Pos = "Arterial_Endo_RAPGEF3+SELP+",
     Arterial_Neg = "Arterial_RAPGEF3-_Cells",
     Venous_Pos = "Venous_Endo_RAPGEF3+SELP-",
     Venous_Neg = "Venous_RAPGEF3-_Cells"
 )
# # Define TF list
my_tfs <- list(
    inflammation=c("NFKB1", "RELA", "STAT1", "STAT3", "IRF1", "IRF7", "JUN", "FOS"),
    vasoconstriction=c("GATA2", "GATA3", "SRF", "MEF2C", "ETV1"),
    cAMP=c("CREB1", "ATF1", "ATF4")
)
#
# # Run the function with optional modules enabled and mistyR per domain
# saved_rds_path <- run_spatial_analysis_modular(
#     data_dir = spatial_data_folder,
#     ref_path = reference_rds,
#     output_dir = results_directory,
#     sample_name = sample_id,
#     target_gene = gene_of_interest,
#     cell_type_mapping = cell_mapping,
#     run_banksy = TRUE,          # <<< Enable Banksy
#     run_cellchat = TRUE,        # <<< Enable CellChat
#     run_dorothea = TRUE,        # <<< Enable Dorothea
#     mistyR_per_domain = TRUE,   # <<< Enable mistyR per domain
#     dorothea_tfs_list = my_tfs,
#     plot_file_format = "pdf"    # <<< Save plots as PDF
# )
#
# # Example: Run only Seurat + RCTD + Global mistyR (no optional modules)
# # saved_rds_path_basic <- run_spatial_analysis_modular(
# #     data_dir = spatial_data_folder,
# #     ref_path = reference_rds,
# #     output_dir = "Analysis_Results/FX0028_Basic_Modular",
# #     sample_name = "FX0028_Basic",
# #     run_banksy = FALSE,
# #     run_cellchat = FALSE,
# #     run_dorothea = FALSE,
# #     mistyR_per_domain = FALSE,
# #     plot_file_format = "pdf"
# # )
#
# if (!is.null(saved_rds_path) && file.exists(saved_rds_path)) {
#     print(paste("Analysis successful. Final Seurat object saved at:", saved_rds_path))
# } else {
#     print("Analysis encountered errors or final object not saved.")
# }
misty_l_param = 10; misty_cutoff_intra = 0.5; misty_cutoff_para = 0.3;
n_cores = max(1, future::availableCores() - 2); rctd_cores = 4;
sct_assay = "SCT"; rna_assay = "SCT";
ref_ident_col = "updated_celltype"; ref_counts_col = "nCount_RNA";
banksy_k_geom = 15; banksy_lambda = 0.2; banksy_k_centers = 10;
plot_file_format = "pdf"

