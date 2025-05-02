# Refactored R Pipeline for cNMF and hdWGCNA Integration
# Description: This script provides a structured workflow for single-cell RNA-seq
# analysis combining cNMF, functional annotation (GSEA), and hdWGCNA.
# It is refactored into functions for better modularity and maintainability.
# Includes option to load a pre-existing Seurat object and specify reduction names.
# Includes an attempt to automatically select optimal K for cNMF.
# Version: Fixes v4 - Added auto K selection attempt.

# --- 0. Load Libraries ---
# Ensure all necessary packages are installed before running.
# install.packages(c("SingleCellExperiment", "Matrix", "reticulate", "Seurat", "SeuratObject",
#                    "hdWGCNA", "tidyverse", "clusterProfiler", "msigdbr",
#                    "enrichplot", "pheatmap", "ComplexHeatmap", "viridis",
#                    "patchwork", "anndata", "BiocParallel", "rlang")) # Add any missing ones

load_libraries <- function() {
  message("Loading required libraries...")
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
    library(reticulate)
    library(Seurat)
    library(SeuratObject) # Explicitly load for class checks if needed
    library(hdWGCNA)
    library(tidyverse) # Includes dplyr, ggplot2, etc.
    library(clusterProfiler)
    library(msigdbr)
    library(enrichplot)
    library(pheatmap)
    library(ComplexHeatmap)
    library(viridis)
    library(patchwork)
    library(tibble)
    library(BiocParallel) # For potentially parallelizing GSEA
    library(rlang) # Needed for %||% operator
  })
  message("Libraries loaded.")
}

# --- 1. Setup and Preprocessing ---
# (Code remains the same as v3)
setup_analysis <- function(config) {
  message("--- 1. Setting up Analysis ---")

  # 1a. Configure Python Environment (Optional)
  if (!is.null(config$python_env_path)) {
    message("Configuring reticulate Python environment: ", config$python_env_path)
    tryCatch({
      if (config$python_env_type == "conda") {
        reticulate::use_condaenv(config$python_env_path, required = TRUE)
      } else if (config$python_env_type == "virtualenv") {
        reticulate::use_virtualenv(config$python_env_path, required = TRUE)
      } else {
         reticulate::use_python(config$python_env_path, required = TRUE)
      }
    }, error = function(e) {
      warning("Could not configure Python environment or find cNMF via reticulate. ",
              "Manual cNMF execution will be required. Error: ", e$message)
    })
  } else {
     message("Python environment path not provided. Assuming manual cNMF execution.")
  }


  # 1b. Load Data (Conditional based on config$input_is_seurat_object)
  message("Loading input data from: ", config$input_data_path)
  if (!file.exists(config$input_data_path)) {
    stop("Input data file not found: ", config$input_data_path)
  }

  if (config$input_is_seurat_object) {
    # Load directly as Seurat object
    message("Loading pre-existing Seurat object...")
    seurat_obj <- readRDS(config$input_data_path)
    message("Input Seurat object loaded.")

    # --- Checks for pre-loaded Seurat object ---
    if (!inherits(seurat_obj, "Seurat")) {
        stop("Loaded object from ", config$input_data_path, " is not a Seurat object.")
    }
    message("Object class verified as Seurat.")

    # Check for cell type annotation
    if (!config$celltype_col %in% colnames(seurat_obj@meta.data)) {
        stop("Cell type column '", config$celltype_col, "' not found in the metadata of the loaded Seurat object.")
    }
    message("Cell type column '", config$celltype_col, "' found.")

    # Check for raw counts assay (needed for cNMF)
    default_assay <- DefaultAssay(seurat_obj)
    if (!"counts" %in% Layers(seurat_obj[[default_assay]])) {
         found_counts <- FALSE
         assay_list <- names(seurat_obj@assays)
         if (!is.null(assay_list) && length(assay_list) > 0) {
             for (assay_name in assay_list) {
                 if ("counts" %in% Layers(seurat_obj[[assay_name]])) {
                     counts_data <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
                     if (!is.null(counts_data) && nrow(counts_data) > 0 && ncol(counts_data) > 0) {
                        message("Found 'counts' layer in assay: ", assay_name, ". This will be used for cNMF.")
                        found_counts <- TRUE
                        break
                     }
                 }
             }
         } else {
             warning("Could not retrieve assay names using names(seurat_obj@assays). Cannot search for 'counts' layer.")
         }

         if (!found_counts) {
            warning("Could not find a 'counts' layer with data in any assay of the loaded Seurat object. ",
                    "The 'prepare_cnmf_input' step might fail as cNMF requires raw counts.")
         }
    } else {
         message("Raw 'counts' layer found in default assay '", default_assay, "'.")
    }


    # Check if user-specified reductions exist
    pca_reduction_name <- config$seurat_input_pca_name %||% "pca" # Default to "pca" if NULL
    umap_reduction_name <- config$seurat_input_umap_name %||% "umap" # Default to "umap" if NULL

    if (!pca_reduction_name %in% names(seurat_obj@reductions)) {
        warning("PCA reduction named '", pca_reduction_name, "' not found in the loaded Seurat object. Some downstream functions might expect it.")
    } else {
        message("PCA reduction named '", pca_reduction_name, "' found.")
    }
     if (!umap_reduction_name %in% names(seurat_obj@reductions)) {
        warning("UMAP reduction named '", umap_reduction_name, "' not found in the loaded Seurat object. Some downstream functions might expect it.")
    } else {
        message("UMAP reduction named '", umap_reduction_name, "' found.")
    }

    message("Skipping standard Seurat preprocessing steps as requested.")

  } else {
    # Load as SingleCellExperiment and preprocess
    message("Loading as SingleCellExperiment object and performing preprocessing...")
    sce <- readRDS(config$input_data_path)
    message("Input SCE data loaded.")

    if (!inherits(sce, "SingleCellExperiment")) {
        stop("Loaded object from ", config$input_data_path, " is not a SingleCellExperiment object. Set 'input_is_seurat_object = TRUE' if it's a Seurat object.")
    }

    # Check for raw counts in SCE
    if (!"counts" %in% names(assays(sce))) {
        stop("Raw counts are required in the 'counts' assay slot of the sce object.")
    }
    # Check for cell type annotation in SCE
    if (!config$celltype_col %in% colnames(colData(sce))) {
        stop("Cell type column '", config$celltype_col, "' not found in colData(sce).")
    }

    # Create and Preprocess Seurat Object
    message("Creating Seurat object from SCE...")
    seurat_obj <- Seurat::CreateSeuratObject(counts = assays(sce)$counts,
                                             meta.data = as.data.frame(colData(sce)))

    message("Preprocessing Seurat object (Normalize, Find HVGs, Scale, PCA, UMAP)...")
    seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = FALSE, nfeatures = config$seurat_n_hvgs)
    seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE, npcs = config$seurat_n_pcs)
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:config$seurat_umap_dims, verbose = FALSE)
    message("Seurat object created and preprocessed.")
    # Store the standard reduction names used during preprocessing
    pca_reduction_name <- "pca"
    umap_reduction_name <- "umap"
  }

  # Final checks and print info
  message("--- Setup Summary ---")
  print(seurat_obj)
  # Ensure DimPlot uses the correct UMAP reduction name and check existence
  # Use umap_reduction_name determined above
  if (exists("umap_reduction_name") && umap_reduction_name %in% names(seurat_obj@reductions)) {
      print(DimPlot(seurat_obj, group.by = config$celltype_col, label = TRUE, reduction = umap_reduction_name) + NoLegend() + ggtitle(paste("UMAP by Cell Type (Reduction:", umap_reduction_name, ")")))
  } else {
      # Try default "umap" if the determined one doesn't exist (e.g., if loaded object had different name but user didn't specify)
      if ("umap" %in% names(seurat_obj@reductions)) {
           print(DimPlot(seurat_obj, group.by = config$celltype_col, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("UMAP by Cell Type (Reduction: umap)"))
      } else {
           message("Default UMAP reduction ('umap') not available for plotting.")
      }
  }
  message("--- End Setup ---")

  return(seurat_obj)
}


# --- 2. cNMF Analysis ---

# 2a. Prepare Input for cNMF (Generates files and instructions)
# (Code remains the same as v3)
prepare_cnmf_input <- function(seurat_obj, config) {
  message("--- 2a. Preparing cNMF Input ---")
  output_dir <- file.path(config$output_basedir, "cnmf_analysis")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  message("cNMF output will be stored in: ", normalizePath(output_dir))

  # --- Get Raw Counts ---
  counts_matrix <- NULL
  counts_assay_name <- NULL
  default_assay <- DefaultAssay(seurat_obj)

  if ("counts" %in% Layers(seurat_obj[[default_assay]])) {
      counts_matrix <- GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
      if (!is.null(counts_matrix) && nrow(counts_matrix) > 0 && ncol(counts_matrix) > 0) {
         counts_assay_name <- default_assay
      } else {
         counts_matrix <- NULL
         message("Found 'counts' layer in default assay '", default_assay,"', but it appears empty.")
      }
  }

  if (is.null(counts_matrix)) {
      message("Searching other assays for non-empty 'counts' layer...")
      assay_list <- names(seurat_obj@assays)
      if (!is.null(assay_list) && length(assay_list) > 0) {
          for (assay_name in assay_list) {
              if (assay_name == default_assay) next
              if ("counts" %in% Layers(seurat_obj[[assay_name]])) {
                  counts_matrix_candidate <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
                  if (!is.null(counts_matrix_candidate) && nrow(counts_matrix_candidate) > 0 && ncol(counts_matrix_candidate) > 0) {
                      counts_matrix <- counts_matrix_candidate
                      counts_assay_name <- assay_name
                      message("Using 'counts' layer from assay: ", counts_assay_name)
                      break
                  } else {
                      message("Found 'counts' layer in assay '", assay_name,"', but it appears empty.")
                  }
              }
          }
      } else {
          warning("Could not retrieve assay names using names(seurat_obj@assays) during counts search in prepare_cnmf_input.")
      }
  }


  if (is.null(counts_matrix)) {
      stop("Could not find a usable 'counts' layer in any assay of the Seurat object. cNMF requires raw counts.")
  }
  # --- End Get Raw Counts ---


  counts_matrix_sparse <- as(counts_matrix, "dgCMatrix")

  # Save as h5ad (recommended)
  h5ad_path <- file.path(output_dir, "counts.h5ad")
  message("Saving counts from assay '", counts_assay_name, "' to .h5ad format: ", h5ad_path)
  tryCatch({
      anndata_py <- reticulate::import("anndata")
      scipy_sparse <- reticulate::import("scipy.sparse")

      # AnnData needs cells x genes
      adata <- anndata_py$AnnData(
          X = scipy_sparse$csc_matrix(t(counts_matrix_sparse)),
          obs = seurat_obj@meta.data,
          var = data.frame(gene_name = rownames(counts_matrix_sparse), row.names = rownames(counts_matrix_sparse))
      )
      adata$write_h5ad(h5ad_path)
      message("Counts saved successfully to counts.h5ad")
  }, error = function(e) {
      warning("Failed to save as .h5ad using reticulate. Check 'anndata' and 'scipy' Python packages. ",
              "Consider saving manually or using a different format if needed. Error: ", e$message)
  })

  # Generate cNMF commands
  k_range_str <- paste(config$cnmf_k_range, collapse = " ")
  prepare_cmd <- paste0("cnmf prepare --output-dir . --name ", config$cnmf_run_name,
                        " -c ", basename(h5ad_path), " -k ", k_range_str,
                        " --n-iter ", config$cnmf_n_iter, " -j ", config$cnmf_num_workers,
                        " --seed ", config$cnmf_seed)

  factorize_cmd_base <- paste0("cnmf factorize --output-dir . --name ", config$cnmf_run_name)
  factorize_cmds <- paste0(factorize_cmd_base, " --worker-index ",
                           0:(config$cnmf_num_workers - 1), " --total-workers ",
                           config$cnmf_num_workers, collapse = "\n")

  combine_cmd <- paste0("cnmf combine --output-dir . --name ", config$cnmf_run_name)

  message("\n--- cNMF EXTERNAL EXECUTION INSTRUCTIONS ---")
  message("cNMF is computationally intensive. Run the following commands in your terminal,")
  message("after navigating to the cNMF output directory:")
  message("cd ", normalizePath(output_dir))
  message("\n1. Prepare:")
  message(prepare_cmd)
  message("\n2. Factorize (run one command per worker, potentially in parallel):")
  message(factorize_cmds)
  message("\n3. Combine:")
  message(combine_cmd)
  message("\n4. K Selection:")
  message("After completion, examine the '", config$cnmf_run_name,
          ".k_selection.png' plot in the output directory to choose the optimal K.")
  message("Set 'cnmf_optimal_k' in the config accordingly OR leave it NULL to attempt automatic selection.")
  message("--- END cNMF INSTRUCTIONS ---\n")


  return(list(output_dir = output_dir, run_name = config$cnmf_run_name))
}

# NEW Function: Attempt to automatically select K based on stability file
select_optimal_k_auto <- function(cnmf_output_dir, cnmf_run_name, k_range) {
  message("--- Attempting Automatic K Selection ---")
  stats_file <- file.path(cnmf_output_dir, paste0(cnmf_run_name, ".stats.df.tsv"))

  if (!file.exists(stats_file)) {
    warning("cNMF stats file not found: ", stats_file, ". Cannot perform automatic K selection.")
    return(NULL)
  }

  tryCatch({
    stats_df <- read.delim(stats_file, sep = "\t")

    # Basic check for expected columns
    if (!all(c("k", "stability") %in% colnames(stats_df))) {
      warning("Stats file '", stats_file, "' does not contain expected 'k' and 'stability' columns.")
      return(NULL)
    }

    # Filter for the tested K range (in case the file contains others)
    stats_df <- stats_df %>% filter(k %in% k_range)

    if (nrow(stats_df) == 0) {
        warning("No stability data found for the specified K range in the stats file.")
        return(NULL)
    }

    # Simple Strategy: Choose K with maximum stability
    # More complex strategies (e.g., elbow method on error) could be added here
    best_k <- stats_df %>%
      arrange(desc(stability)) %>%
      slice(1) %>%
      pull(k)

    message("Automatic K selection based on maximum stability suggests K = ", best_k)
    # Add a check: sometimes max stability is at the lowest K, which might not be ideal.
    # A more robust method might find the peak before a major drop, or use an elbow method.
    # For now, we'll use the max stability K.
    return(as.integer(best_k))

  }, error = function(e) {
    warning("Error reading or processing cNMF stats file '", stats_file, "': ", e$message)
    return(NULL)
  })
}


# 2b. Load cNMF Results
# (Code remains the same as v3)
load_cnmf_results <- function(seurat_obj, cnmf_output_dir, cnmf_run_name, optimal_k) {
  message("--- 2b. Loading cNMF Results for K = ", optimal_k, " ---")

  usage_file <- file.path(cnmf_output_dir, paste0(cnmf_run_name, ".usages.k_", optimal_k, ".consensus.txt"))
  spectra_file <- file.path(cnmf_output_dir, paste0(cnmf_run_name, ".spectra.k_", optimal_k, ".consensus.txt"))

  if (!file.exists(usage_file) || !file.exists(spectra_file)) {
    stop("cNMF output files not found for k=", optimal_k, ". Check the path: ", cnmf_output_dir,
         ", run name: ", cnmf_run_name, " and selected K.")
  }

  # Load H matrix (Usage): Programs x Cells -> transpose to Cells x Programs
  H_matrix_raw <- read.table(usage_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  H_matrix <- t(as.matrix(H_matrix_raw)) # Now Cells x Programs
  colnames(H_matrix) <- paste0("NMF_Program_", 1:ncol(H_matrix))

  # Load W matrix (Spectra): Genes x Programs (assuming this format)
  W_matrix_raw <- read.table(spectra_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  W_matrix <- as.matrix(W_matrix_raw)
  colnames(W_matrix) <- paste0("NMF_Program_", 1:ncol(W_matrix))

  # Verify dimensions and gene names
  # Get features from a relevant assay (e.g., RNA or default)
  target_assay_features <- rownames(seurat_obj) # Default to all features in object
  if ("RNA" %in% names(seurat_obj@assays)) {
      target_assay_features <- rownames(seurat_obj[["RNA"]])
  } else {
      target_assay_features <- rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
  }


  if (nrow(W_matrix) != length(target_assay_features)) {
      warning("Mismatch in number of genes between W matrix (", nrow(W_matrix), ") and Seurat object features (", length(target_assay_features), "). Trying to subset W matrix.")
       common_genes <- intersect(rownames(W_matrix), target_assay_features)
       if (length(common_genes) == 0) stop("No common genes found between W matrix and Seurat object.")
       message("Found ", length(common_genes), " common genes.")
       W_matrix <- W_matrix[common_genes, , drop = FALSE]
       # Reorder W matrix to match Seurat object gene order (among common genes)
       W_matrix <- W_matrix[intersect(target_assay_features, rownames(W_matrix)), ]

  } else if (!all(rownames(W_matrix) %in% target_assay_features)) {
       warning("Some genes in W matrix are not present in the Seurat object features. Subsetting W matrix.")
       common_genes <- intersect(rownames(W_matrix), target_assay_features)
       if (length(common_genes) == 0) stop("No common genes found between W matrix and Seurat object.")
       message("Found ", length(common_genes), " common genes.")
       W_matrix <- W_matrix[common_genes, , drop = FALSE]
       W_matrix <- W_matrix[intersect(target_assay_features, rownames(W_matrix)), ]
   } else {
        # Ensure W matrix gene order matches Seurat object features
        W_matrix <- W_matrix[target_assay_features, ]
   }


  # Verify cell names and order for H matrix
  if (!all(rownames(H_matrix) %in% colnames(seurat_obj))) {
      stop("Some cells in H matrix are not present in the Seurat object.")
  }
  # Ensure H matrix cell order matches Seurat object
  H_matrix <- H_matrix[colnames(seurat_obj), , drop = FALSE]

  # Add NMF program activities (H matrix) to Seurat metadata
  message("Adding NMF program activities to Seurat metadata.")
  seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(H_matrix))

  message("cNMF results loaded and added to Seurat object.")
  return(list(seurat_obj = seurat_obj, W_matrix = W_matrix, H_matrix = H_matrix))
}


# 2c. Visualize cNMF Results
# (Code remains the same as v3)
visualize_cnmf_results <- function(seurat_obj, H_matrix, W_matrix, config) {
  message("--- 2c. Visualizing cNMF Results ---")
  plot_dir <- file.path(config$output_basedir, "plots", "cnmf")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  # UMAPs colored by program scores
  message("Generating UMAP plots for NMF programs...")
  nmf_programs <- colnames(H_matrix)
  n_plots <- min(length(nmf_programs), config$cnmf_num_umap_plots) # Limit plots
  plot_list_umap <- list()
  # Determine UMAP reduction name to use
  umap_reduction_name <- config$seurat_input_umap_name %||% "umap"

  if (umap_reduction_name %in% names(seurat_obj@reductions)) {
      for (program in nmf_programs[1:n_plots]) {
        p <- FeaturePlot(seurat_obj, features = program, coord.fixed = TRUE, pt.size = 0.1, reduction = umap_reduction_name) +
             scale_colour_viridis(option = "plasma") +
             ggtitle(program)
        plot_list_umap[[program]] <- p
      }
      if (length(plot_list_umap) > 0) {
          combined_plot_umap <- wrap_plots(plot_list_umap, ncol = min(3, n_plots)) +
                                plot_annotation(title = paste("NMF Program Activity on UMAP (Reduction:", umap_reduction_name, ")"))
          print(combined_plot_umap)
          ggsave(file.path(plot_dir, "cnmf_umap_program_activity.png"), combined_plot_umap, width = 12, height = 4 * ceiling(n_plots/3), dpi = 300)
      }
  } else {
      message("UMAP reduction '", umap_reduction_name, "' not found, skipping NMF program UMAP plots.")
  }


  # Heatmap of program activities across cells (sampled)
  message("Generating heatmap of NMF program activities...")
  cells_to_plot <- sample(colnames(seurat_obj), min(config$cnmf_heatmap_ncells, ncol(seurat_obj)))
  # Ensure celltype column exists before ordering/subsetting
  if (!config$celltype_col %in% colnames(seurat_obj@meta.data)) {
       message("Celltype column '", config$celltype_col, "' not found. Heatmap will not be ordered by cell type.")
       cell_order_subset <- cells_to_plot
       annotation_col <- NULL
       anno_colors <- NULL
       col_split <- NULL
  } else {
      cell_order <- colnames(seurat_obj)[order(seurat_obj@meta.data[[config$celltype_col]])]
      cell_order_subset <- intersect(cell_order, cells_to_plot)
      meta_subset <- seurat_obj@meta.data[cell_order_subset, , drop = FALSE]
      H_subset <- H_matrix[cell_order_subset, , drop = FALSE]

      annotation_col = data.frame(
        Celltype = meta_subset[[config$celltype_col]]
      )
      rownames(annotation_col) = rownames(meta_subset)

      # Create colors for cell types
      celltypes <- unique(seurat_obj@meta.data[[config$celltype_col]])
      celltype_levels <- levels(factor(seurat_obj@meta.data[[config$celltype_col]]))
      celltype_colors <- scales::hue_pal()(length(celltype_levels))
      names(celltype_colors) <- celltype_levels
      anno_colors <- list(Celltype = celltype_colors)
      col_split <- annotation_col$Celltype # Split columns by cell type
  }

  # Adjust H_subset if annotation_col is NULL
  if(is.null(annotation_col)){
      H_subset <- H_matrix[cell_order_subset, , drop = FALSE]
      heatmap_annotation <- NULL
  } else {
      H_subset <- H_matrix[cell_order_subset, , drop = FALSE]
      heatmap_annotation <- HeatmapAnnotation(df = annotation_col, col = anno_colors)
  }


  # Use ComplexHeatmap for better control
  h_heatmap <- ComplexHeatmap::Heatmap(t(H_subset), # Programs x Cells
                      name = "Scaled Activity",
                      column_title = paste("NMF Program Activity (Sampled", length(cell_order_subset), "Cells)"),
                      row_title = "NMF Programs",
                      show_column_names = FALSE,
                      cluster_columns = is.null(annotation_col), # Cluster cols only if no annotation/split
                      cluster_rows = TRUE,
                      col = viridis(100),
                      heatmap_legend_param = list(title = "Scaled Activity"),
                      top_annotation = heatmap_annotation,
                      row_names_gp = gpar(fontsize = 8),
                      column_split = col_split
                      )

  png(file.path(plot_dir, "cnmf_heatmap_program_activity.png"), width = 10, height = 8, units = "in", res = 300)
  draw(h_heatmap)
  dev.off()
  print("Heatmap saved.")


  # Heatmap of top genes per program (W matrix)
  # (Code remains the same as v3)
  message("Generating heatmap of top genes per NMF program...")
  n_top_genes <- config$cnmf_num_top_genes_heatmap
  top_genes_list <- list()
  W_matrix_df <- as.data.frame(W_matrix) # Use potentially subsetted W_matrix

  for (program in colnames(W_matrix_df)) {
      if (!program %in% names(W_matrix_df)) next
      top_genes_data <- W_matrix_df %>%
          rownames_to_column("gene") %>%
          dplyr::select(gene, !!sym(program)) %>%
          arrange(desc(!!sym(program)))
      top_n_values <- head(unique(top_genes_data[[program]]), n_top_genes)
      min_cutoff <- if(length(top_n_values) > 0) min(top_n_values) else -Inf
      top_genes <- top_genes_data %>%
                   filter(!!sym(program) >= min_cutoff) %>%
                   arrange(desc(!!sym(program))) %>%
                   slice_head(n = n_top_genes) %>%
                   pull(gene)
      top_genes_list[[program]] <- top_genes
  }

  genes_to_plot_w <- unique(unlist(top_genes_list))
  if (length(genes_to_plot_w) > 0) {
      genes_to_plot_w <- intersect(genes_to_plot_w, rownames(W_matrix))
      if (length(genes_to_plot_w) > 0) {
          W_subset_top <- W_matrix[genes_to_plot_w, , drop = FALSE]
          w_heatmap <- ComplexHeatmap::Heatmap(W_subset_top,
                              name = "Gene Loading",
                              show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
                              column_title = paste("Top", n_top_genes, "Genes per NMF Program (W Matrix)"),
                              column_names_rot = 45,
                              cluster_columns = TRUE,
                              cluster_rows = TRUE,
                              row_dend_reorder = TRUE,
                              column_dend_reorder = TRUE)
          png(file.path(plot_dir, "cnmf_heatmap_top_genes.png"), width = 10, height = max(8, length(genes_to_plot_w)*0.1), units = "in", res = 300)
          draw(w_heatmap)
          dev.off()
          print("Top genes heatmap saved.")
      } else {
           message("No top genes found in the current W matrix to plot.")
      }
  } else {
      message("No top genes found to plot for W matrix heatmap.")
  }

  message("cNMF visualization complete.")
}


# --- 3. Functional Annotation (GSEA) ---
# (Code remains the same as v3)
# 3a. Run GSEA on NMF Programs
run_gsea_on_nmf <- function(W_matrix, config) {
  message("--- 3a. Running GSEA for NMF Programs ---")

  # Prepare gene sets
  message("Fetching gene sets from msigdbr for species: ", config$gsea_species)
  gsea_categories <- strsplit(config$gsea_categories, ",")[[1]] # e.g., "H,C5"
  gsea_subcategories <- if(!is.null(config$gsea_subcategories)) strsplit(config$gsea_subcategories, ",")[[1]] else NULL # e.g., "GO:BP"

  msigdbr_collections <- msigdbr(species = config$gsea_species)

  term2gene_list <- list()
  if ("H" %in% gsea_categories) {
     term2gene_list$Hallmark <- msigdbr_collections %>%
        filter(gs_cat == "H") %>%
        dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
     message("Loaded Hallmark gene sets.")
  }
   if ("C5" %in% gsea_categories) {
      if (!is.null(gsea_subcategories) && "GO:BP" %in% gsea_subcategories) {
          term2gene_list$GO_BP <- msigdbr_collections %>%
             filter(gs_cat == "C5", gs_subcat == "GO:BP") %>%
             dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
          message("Loaded GO:BP gene sets.")
      }
      if (!is.null(gsea_subcategories) && "GO:MF" %in% gsea_subcategories) {
           term2gene_list$GO_MF <- msigdbr_collections %>%
             filter(gs_cat == "C5", gs_subcat == "GO:MF") %>%
             dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
           message("Loaded GO:MF gene sets.")
       }
       if (!is.null(gsea_subcategories) && "GO:CC" %in% gsea_subcategories) {
           term2gene_list$GO_CC <- msigdbr_collections %>%
             filter(gs_cat == "C5", gs_subcat == "GO:CC") %>%
             dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
           message("Loaded GO:CC gene sets.")
       }
   }
    if ("C2" %in% gsea_categories) {
        if (!is.null(gsea_subcategories) && "CP:KEGG" %in% gsea_subcategories) {
            term2gene_list$KEGG <- msigdbr_collections %>%
                filter(gs_cat == "C2", gs_subcat == "CP:KEGG") %>%
                dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
            message("Loaded KEGG gene sets.")
        }
        if (!is.null(gsea_subcategories) && "CP:REACTOME" %in% gsea_subcategories) {
            term2gene_list$Reactome <- msigdbr_collections %>%
                filter(gs_cat == "C2", gs_subcat == "CP:REACTOME") %>%
                dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
            message("Loaded Reactome gene sets.")
        }
    }
    if ("C7" %in% gsea_categories) {
         if (!is.null(gsea_subcategories) && "IMMUNESIGDB" %in% gsea_subcategories) {
            term2gene_list$ImmuneSigDB <- msigdbr_collections %>%
                filter(gs_cat == "C7", gs_subcat == "IMMUNESIGDB") %>%
                dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
            message("Loaded ImmuneSigDB gene sets.")
        }
    }


   if (length(term2gene_list) == 0) {
       stop("No valid gene set categories specified or loaded.")
   }

  # Run GSEA for each program and each gene set collection
  gsea_results_list <- list()
  programs <- colnames(W_matrix) # Use potentially subsetted W_matrix

  for (program in programs) {
    message("Running GSEA for ", program)
    gene_list <- W_matrix[, program]
    gene_list <- sort(gene_list, decreasing = TRUE)
    gene_list <- gene_list[gene_list != 0] # Remove zero scores

    if(length(gene_list) == 0) {
        message("  Skipping GSEA for ", program, ": No non-zero gene scores.")
        gsea_results_list[[program]] <- list()
        next
    }
     if(is.null(names(gene_list))) {
        message("  Skipping GSEA for ", program, ": Gene list is missing names.")
        gsea_results_list[[program]] <- list()
        next
    }
    if(any(duplicated(names(gene_list)))) {
        message("  Warning: Duplicate gene names found for program ", program, ". Keeping highest score.")
        gene_list <- gene_list[!duplicated(names(gene_list))]
    }


    program_gsea_results <- list()
    for (gs_name in names(term2gene_list)) {
       message("  Gene Set: ", gs_name)
       term2gene <- term2gene_list[[gs_name]]
       tryCatch({
         gsea_res <- GSEA(geneList = gene_list,
                          TERM2GENE = term2gene,
                          pvalueCutoff = config$gsea_pvalue_cutoff,
                          pAdjustMethod = "BH",
                          verbose = FALSE,
                          minGSSize = config$gsea_min_geneset_size,
                          maxGSSize = config$gsea_max_geneset_size,
                          seed = TRUE,
                          eps = 0)

         if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
            program_gsea_results[[gs_name]] <- as.data.frame(gsea_res@result)
         } else {
            program_gsea_results[[gs_name]] <- data.frame()
         }

       }, error = function(e) {
         message("   GSEA failed for ", program, " with gene set ", gs_name, ": ", e$message)
         program_gsea_results[[gs_name]] <- data.frame()
       })
    }
     gsea_results_list[[program]] <- program_gsea_results
  }

  message("GSEA analysis complete.")
  return(gsea_results_list)
}

# 3b. Visualize GSEA Results
# (Code remains the same as v3)
visualize_gsea_results <- function(gsea_results_list, config) {
  message("--- 3b. Visualizing GSEA Results ---")
  plot_dir <- file.path(config$output_basedir, "plots", "gsea")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  if (length(gsea_results_list) == 0 || all(sapply(gsea_results_list, length) == 0)) {
      message("No GSEA results found to visualize.")
      return()
  }

  first_valid_result_index <- which(sapply(gsea_results_list, function(x) !is.null(x) && length(x) > 0))[1]
  if (is.na(first_valid_result_index)) {
       message("Could not find any valid GSEA results to determine gene set names.")
       return()
  }
  first_valid_result <- gsea_results_list[[first_valid_result_index]]
  if (is.null(first_valid_result) || length(first_valid_result) == 0) {
      message("Could not determine gene set names from GSEA results.")
      return()
  }
  gene_set_names <- names(first_valid_result)

  for (gs_name in gene_set_names) {
    message("Generating GSEA summary plot for: ", gs_name)
    all_gsea_results_for_gs <- bind_rows(
      lapply(names(gsea_results_list), function(prog_name) {
        if (!is.null(gsea_results_list[[prog_name]]) && gs_name %in% names(gsea_results_list[[prog_name]])) {
            res_df <- gsea_results_list[[prog_name]][[gs_name]]
            if (!is.null(res_df) && nrow(res_df) > 0) {
              res_df %>%
                mutate(Program = prog_name) %>%
                filter(p.adjust < config$gsea_plot_fdr_cutoff)
            } else { NULL }
        } else { NULL }
      })
    )

    if (!is.null(all_gsea_results_for_gs) && nrow(all_gsea_results_for_gs) > 0) {
        top_terms_per_program <- all_gsea_results_for_gs %>%
             mutate(Description_clean = sub("^(HALLMARK|GOBP|GOMF|GOCC|KEGG|REACTOME)_", "", Description)) %>%
            group_by(Program) %>%
            arrange(p.adjust) %>%
            slice_head(n = config$gsea_plot_top_n_terms) %>%
            ungroup() %>%
             mutate(Description_clean = factor(Description_clean, levels = unique(Description_clean[order(NES, decreasing = TRUE)])))

        gsea_dotplot <- ggplot(top_terms_per_program,
                               aes(x = Program, y = Description_clean, size = setSize, color = NES)) +
            geom_point() +
            scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name="NES") +
            scale_size(range = c(2, 8), name="Set Size") +
            theme_bw(base_size = 10) +
            ylab(paste(gs_name, "Gene Sets")) +
            xlab("NMF Program") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  axis.text.y = element_text(size=8)) +
            ggtitle(paste("Top", config$gsea_plot_top_n_terms, "Significant", gs_name, "GSEA Results per NMF Program"))

        print(gsea_dotplot)
        ggsave(file.path(plot_dir, paste0("gsea_summary_dotplot_", gs_name, ".png")),
               gsea_dotplot, width = max(10, length(unique(top_terms_per_program$Program))*0.8),
               height = max(6, length(unique(top_terms_per_program$Description_clean))*0.3), dpi = 300, limitsize = FALSE)

    } else {
        message("No significant GSEA results found for ", gs_name, " at p.adjust < ", config$gsea_plot_fdr_cutoff)
    }
  }
  message("GSEA visualization complete.")
}


# --- 4. hdWGCNA Analysis ---
# 4a. Setup hdWGCNA Environment
# (Code remains the same as v3)
setup_hdwgcna <- function(seurat_obj, config) {
  message("--- 4a. Setting up hdWGCNA ---")

  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = config$hdwgcna_gene_select_method,
    features = if(config$hdwgcna_gene_select_method == "custom") config$hdwgcna_custom_features else NULL,
    wgcna_name = config$hdwgcna_name
  )

  if (config$hdwgcna_use_metacells) {
    message("Constructing metacells by group: ", paste(config$hdwgcna_metacell_group_by, collapse=", "))
    if (!all(config$hdwgcna_metacell_group_by %in% colnames(seurat_obj@meta.data))) {
        stop("Metacell grouping variable(s) not found in Seurat metadata: ",
             paste(config$hdwgcna_metacell_group_by[!config$hdwgcna_metacell_group_by %in% colnames(seurat_obj@meta.data)], collapse=", "))
    }
    seurat_obj <- MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = config$hdwgcna_metacell_group_by,
      k = config$hdwgcna_metacell_k,
      max_shared = config$hdwgcna_metacell_max_shared,
      ident.group = config$hdwgcna_metacell_ident
    )
    metacell_assay <- paste0(config$hdwgcna_assay_for_wgcna, "_metacell")
    if (!metacell_assay %in% names(seurat_obj@assays)) {
        stop("Metacell assay '", metacell_assay, "' not found after running MetacellsByGroups.")
    }
    seurat_obj <- NormalizeData(seurat_obj, assay = metacell_assay, normalization.method = 'LogNormalize', verbose=FALSE)
    seurat_obj <- ScaleData(seurat_obj, assay = metacell_assay, verbose=FALSE)

    message("Setting expression data for WGCNA using metacells.")
    seurat_obj <- SetDatExpr(
      seurat_obj,
      group_name = config$hdwgcna_group_name,
      group.by = config$hdwgcna_metacell_group_by,
      assay = metacell_assay,
      layer = config$hdwgcna_layer_for_wgcna
    )
     seurat_obj@misc[[config$hdwgcna_name]]$metacell_ident <- config$hdwgcna_metacell_ident
  } else {
    message("Setting expression data for WGCNA using single cells.")
     if (!is.null(config$hdwgcna_group.by) && !all(config$hdwgcna_group.by %in% colnames(seurat_obj@meta.data))) {
        stop("hdWGCNA grouping variable(s) not found in Seurat metadata: ",
             paste(config$hdwgcna_group.by[!config$hdwgcna_group.by %in% colnames(seurat_obj@meta.data)], collapse=", "))
    }
    seurat_obj <- SetDatExpr(
      seurat_obj,
      group_name = config$hdwgcna_group_name,
      group.by = config$hdwgcna_group.by,
      assay = config$hdwgcna_assay_for_wgcna,
      layer = config$hdwgcna_layer_for_wgcna
    )
    seurat_obj@misc[[config$hdwgcna_name]]$metacell_ident <- NULL
  }

  message("Testing soft power thresholds...")
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = config$hdwgcna_network_type,
    setDatExpr = FALSE
  )

  plot_dir <- file.path(config$output_basedir, "plots", "hdwgcna")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  message("Plotting soft power thresholds...")
  plot_list <- PlotSoftPowers(seurat_obj)
  softpower_plot <- wrap_plots(plot_list, ncol = 2) + plot_annotation(title = "hdWGCNA Soft Power Selection")
  print(softpower_plot)
  ggsave(file.path(plot_dir, "hdwgcna_softpower_selection.png"), softpower_plot, width = 10, height = 5, dpi = 300)

  power_table <- GetPowerTable(seurat_obj)
  message("Soft power table:")
  print(power_table)
  message("Examine the plot and table to choose the appropriate 'hdwgcna_soft_power' value in the config.")
  message("Typically choose the lowest power where the scale-free topology fit (SFT.R.sq) plateaus above ", config$hdwgcna_min_rsquared)

  message("hdWGCNA setup complete.")
  return(seurat_obj)
}

# 4b. Construct Network and Find Modules
# (Code remains the same as v3)
run_hdwgcna_network <- function(seurat_obj, config) {
  message("--- 4b. Running hdWGCNA Network Construction ---")
  if (is.null(config$hdwgcna_soft_power)) {
      stop("Please set 'hdwgcna_soft_power' in the configuration based on the soft power plots.")
  }
  message("Using soft power: ", config$hdwgcna_soft_power)

  message("Constructing co-expression network (TOM)...")
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power = config$hdwgcna_soft_power,
    setDatExpr = FALSE,
    tom_name = config$hdwgcna_tom_name,
    overwrite_tom = TRUE,
    networkType = config$hdwgcna_network_type
  )

  plot_dir <- file.path(config$output_basedir, "plots", "hdwgcna")
  png(file.path(plot_dir, "hdwgcna_dendrogram.png"), width = 10, height = 6, units = "in", res = 300)
  PlotDendrogram(seurat_obj, main = paste(config$hdwgcna_name, 'Dendrogram'))
  dev.off()
  print("Dendrogram saved.")

  message("Computing module eigengenes (MEs)...")
  group.by.var <- seurat_obj@misc[[config$hdwgcna_name]]$metacell_ident
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = group.by.var
  )

  message("Computing module connectivity (kME)...")
  kme_group.by <- if (!is.null(group.by.var)) config$hdwgcna_metacell_ident else NULL
  kme_group_name <- if (!is.null(group.by.var)) config$hdwgcna_metacell_ident else config$hdwgcna_group_name

  seurat_obj <- ModuleConnectivity(
      seurat_obj,
      group.by = kme_group.by,
      group_name = kme_group_name
  )

  if (config$hdwgcna_compute_module_scores) {
      message("Computing module expression scores...")
      umap_reduction_name <- config$seurat_input_umap_name %||% "umap"
      if (!umap_reduction_name %in% names(seurat_obj@reductions)) {
          warning("UMAP reduction '", umap_reduction_name, "' not found. Module scores will be computed but cannot be plotted on UMAP by default.")
      }
      seurat_obj <- ModuleExprScore(
        seurat_obj,
        n_genes = config$hdwgcna_module_score_n_genes,
        method = 'Seurat'
      )
  }

  message("hdWGCNA network construction and module definition complete.")
  return(seurat_obj)
}


# 4c. Visualize hdWGCNA Modules
# (Code remains the same as v3)
visualize_hdwgcna_modules <- function(seurat_obj, config) {
  message("--- 4c. Visualizing hdWGCNA Modules ---")
  plot_dir <- file.path(config$output_basedir, "plots", "hdwgcna")
  # Determine UMAP reduction name to use
  umap_reduction_name <- config$seurat_input_umap_name %||% "umap"

  if (config$hdwgcna_compute_module_scores && umap_reduction_name %in% names(seurat_obj@reductions)) {
      message("Generating UMAP plots for module scores...")
      module_colors <- GetModules(seurat_obj)$module_color
      module_colors <- levels(factor(module_colors[module_colors != 'grey']))
      module_names <- paste0("ME", module_colors)

      n_plots <- min(length(module_names), config$hdwgcna_num_umap_plots)
      plot_list_mod_umap <- list()

      score_suffix <- "_Score"
      available_scores <- colnames(seurat_obj@meta.data)

      plotted_count = 0
      for(m in module_names){
          score_name <- paste0(m, score_suffix)
          if (score_name %in% available_scores && plotted_count < n_plots) {
              p <- FeaturePlot(seurat_obj, features = score_name, coord.fixed = TRUE, pt.size = 0.1, reduction = umap_reduction_name) +
                   scale_color_viridis(option="magma") +
                   ggtitle(paste(m, "Score"))
              plot_list_mod_umap[[m]] <- p
              plotted_count <- plotted_count + 1
          } else if (!score_name %in% available_scores) {
              message("  Skipping UMAP for ", m, ": Score '", score_name, "' not found in metadata.")
          }
      }

      if (length(plot_list_mod_umap) > 0) {
          combined_plot_mod_umap <- wrap_plots(plot_list_mod_umap, ncol = min(3, n_plots)) +
                                    plot_annotation(title = paste("hdWGCNA Module Scores on UMAP (Reduction:", umap_reduction_name, ")"))
          print(combined_plot_mod_umap)
          ggsave(file.path(plot_dir, "hdwgcna_umap_module_scores.png"), combined_plot_mod_umap, width = 12, height = 4 * ceiling(plotted_count/3), dpi = 300)
      } else {
          message("No module scores found or plotted.")
      }
  } else if (!config$hdwgcna_compute_module_scores) {
      message("Skipping module score UMAP visualization as 'hdwgcna_compute_module_scores' is FALSE.")
  } else if (!umap_reduction_name %in% names(seurat_obj@reductions)) {
       message("Skipping module score UMAP visualization as UMAP reduction '", umap_reduction_name, "' is not present.")
  }

  message("hdWGCNA module visualization complete.")
}


# --- 5. Integration: NMF Programs and WGCNA Modules ---
# (Code remains the same as v3)
# 5a. Correlate Module Eigengenes with NMF Program Activity
correlate_nmf_hdwgcna <- function(seurat_obj, config) {
  message("--- 5a. Correlating NMF Programs with hdWGCNA Modules ---")

  nmf_trait_names <- grep(paste0("^", config$nmf_program_prefix), colnames(seurat_obj@meta.data), value = TRUE)
  if (length(nmf_trait_names) == 0) {
    stop("NMF program scores (prefixed with '", config$nmf_program_prefix, "') not found in Seurat metadata.")
  }

  metacell_ident <- seurat_obj@misc[[config$hdwgcna_name]]$metacell_ident
  use_harmonized_mes <- !is.null(metacell_ident)

  MEs <- GetMEs(seurat_obj, harmonized = use_harmonized_mes)
  if (is.null(MEs) || nrow(MEs) == 0) {
      stop("Module Eigengenes (MEs) could not be retrieved. Check hdWGCNA steps.")
  }

  if (use_harmonized_mes) {
    message("Aggregating NMF traits to metacell level (", metacell_ident, ")...")
    if (!metacell_ident %in% colnames(seurat_obj@meta.data)) {
        stop("Metacell identifier '", metacell_ident, "' not found in Seurat metadata for aggregation.")
    }
    trait_data_agg <- seurat_obj@meta.data %>%
        rownames_to_column("cell_id") %>%
        dplyr::select(cell_id, !!sym(metacell_ident), all_of(nmf_trait_names)) %>%
        group_by(!!sym(metacell_ident)) %>%
        summarise(across(all_of(nmf_trait_names), mean, na.rm = TRUE)) %>%
        column_to_rownames(metacell_ident)

    if (!all(rownames(MEs) %in% rownames(trait_data_agg))) {
        stop("Some metacell IDs in MEs are missing from aggregated trait data.")
    }
    trait_data <- trait_data_agg[rownames(MEs), , drop = FALSE]
  } else {
    message("Using single-cell NMF traits.")
    if (!all(rownames(MEs) %in% rownames(seurat_obj@meta.data))) {
         stop("Some cell IDs in MEs are missing from Seurat metadata.")
    }
    trait_data <- seurat_obj@meta.data[rownames(MEs), nmf_trait_names, drop = FALSE]
  }

  if(nrow(MEs) != nrow(trait_data)){
      stop("Mismatch in number of observations between MEs (", nrow(MEs), ") and Trait data (", nrow(trait_data), ").")
  }
  if(!all(rownames(MEs) == rownames(trait_data))){
      if(length(intersect(rownames(MEs), rownames(trait_data))) == nrow(MEs)) {
          message("Reordering trait data to match ME observation names.")
          trait_data <- trait_data[rownames(MEs), , drop = FALSE]
      } else {
          stop("Mismatch in observation names (cells/metacells) between MEs and Trait data.")
      }
  }
  message("Dimensions match: ", nrow(MEs), " observations (cells/metacells).")

  message("Computing module-trait correlations...")
  seurat_obj <- ModuleTraitCorrelation(
    seurat_obj,
    traits = trait_data,
    group.by = NULL
  )

  message("Module-trait correlation complete.")
  return(seurat_obj)
}


# 5b. Visualize Module-Trait Correlations
# (Code remains the same as v3)
visualize_module_trait_corr <- function(seurat_obj, config) {
  message("--- 5b. Visualizing Module-Trait Correlations ---")
  plot_dir <- file.path(config$output_basedir, "plots", "integration")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  corr_results <- GetModuleTraitCorrelation(seurat_obj)
  if (is.null(corr_results) || is.null(corr_results$cor) || nrow(corr_results$cor) == 0) {
      message("No module-trait correlation results found to plot.")
      return()
  }

  message("Generating module-trait correlation heatmap...")
  if (is.null(corr_results$cor)) {
       message("Correlation matrix ('cor') is NULL in GetModuleTraitCorrelation results.")
       return()
  }

  corr_plot <- PlotModuleTraitCorrelation(
    seurat_obj,
    label = 'p.value',
    label_symbol = 'stars',
    text_size = 2.5,
    text_digits = 2,
    text_color = 'black',
    high_color = config$integration_plot_high_color,
    mid_color = config$integration_plot_mid_color,
    low_color = config$integration_plot_low_color,
    plot_max = config$integration_plot_corr_max,
    plot_min = -config$integration_plot_corr_max,
    legend_position = 'top',
    legend_title = "Correlation"
  ) + ggtitle("hdWGCNA Module - NMF Program Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          plot.title = element_text(hjust = 0.5))

  print(corr_plot)
  n_mods <- nrow(corr_results$cor)
  n_traits <- ncol(corr_results$cor)
  ggsave(file.path(plot_dir, "hdwgcna_nmf_correlation_heatmap.png"), corr_plot,
         width = max(6, n_traits * 0.5), height = max(5, n_mods * 0.3), dpi = 300, limitsize = FALSE)

  message("Module-trait correlation visualization complete.")
}


# --- 6. Specific Gene Analysis ---
# (Code remains the same as v3)
analyze_gene_nmf_relationship <- function(seurat_obj, W_matrix, config) {
  gene_of_interest <- config$analyze_gene_name
  if (is.null(gene_of_interest) || gene_of_interest == "") {
      message("--- 6. Skipping Specific Gene Analysis (no gene specified) ---")
      return()
  }
  message("--- 6. Analyzing Relationship of Gene: ", gene_of_interest, " with NMF Programs ---")
  plot_dir <- file.path(config$output_basedir, "plots", "gene_analysis", gene_of_interest)
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  gene_found_in_assays <- FALSE
  target_assay_for_expr <- DefaultAssay(seurat_obj)
  assay_list_for_gene_check <- names(seurat_obj@assays)

  if (gene_of_interest %in% rownames(seurat_obj[[target_assay_for_expr]])) {
      gene_found_in_assays <- TRUE
  } else if (!is.null(assay_list_for_gene_check)) {
      for (assay_name in assay_list_for_gene_check) {
          if (gene_of_interest %in% rownames(seurat_obj[[assay_name]])) {
              target_assay_for_expr <- assay_name
              gene_found_in_assays <- TRUE
              message("Gene '", gene_of_interest, "' found in assay: ", target_assay_for_expr)
              break
          }
      }
  }

  if (!gene_found_in_assays) {
      warning("Gene '", gene_of_interest, "' not found in any assay of the Seurat object. Skipping analysis.")
      return()
  }

  if (!gene_of_interest %in% rownames(W_matrix)) {
     warning("Gene '", gene_of_interest, "' not found in the final NMF W matrix (check if it was filtered or in common genes). Skipping analysis.")
     return()
  }

  # 6a. Analyze Gene Contribution (W matrix Loadings)
  message("Analyzing contribution (loading) in W matrix...")
  gene_loadings <- W_matrix[gene_of_interest, ]
  loading_df <- data.frame(
    Program = names(gene_loadings),
    Loading = gene_loadings
  ) %>%
    mutate(Program = factor(Program, levels = names(sort(gene_loadings, decreasing = TRUE))))

  p_loadings <- ggplot(loading_df, aes(x = Program, y = Loading, fill = Loading)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(name = "Loading", direction = -1) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = paste("Contribution (Loading) of", gene_of_interest, "to NMF Programs"),
         x = "NMF Program", y = "Gene Loading (W matrix)")
  print(p_loadings)
  ggsave(file.path(plot_dir, paste0(gene_of_interest, "_nmf_loadings.png")), p_loadings, width = 8, height = 5, dpi = 300)


  # 6b. Analyze Expression Correlation (vs H matrix Activity)
  message("Analyzing correlation with NMF program activity (H matrix)...")
  gene_expression <- GetAssayData(seurat_obj, assay = target_assay_for_expr, layer = "data")[gene_of_interest, ]
  program_colnames <- grep(paste0("^", config$nmf_program_prefix), colnames(seurat_obj@meta.data), value = TRUE)
  program_scores <- seurat_obj@meta.data[, program_colnames, drop = FALSE]

  correlations <- apply(program_scores, 2, function(activity) {
    if (sd(activity, na.rm=TRUE) == 0 || sd(gene_expression, na.rm=TRUE) == 0) return(NA)
    cor(gene_expression, activity, method = "spearman", use="complete.obs")
  })
  correlations[is.na(correlations)] <- 0

  correlation_df <- data.frame(
    Program = names(correlations),
    Correlation = correlations
  ) %>%
    mutate(Program = factor(Program, levels = names(correlations[order(abs(correlations), decreasing = TRUE)])))

  p_correlation <- ggplot(correlation_df, aes(x = Program, y = Correlation, fill = Correlation)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Spearman\nCorr.") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = paste("Correlation:", gene_of_interest, "Expression vs. NMF Program Activity"),
         x = "NMF Program", y = "Spearman Correlation") +
    geom_hline(yintercept = 0, linetype="dashed", color="grey50")
  print(p_correlation)
  ggsave(file.path(plot_dir, paste0(gene_of_interest, "_nmf_correlation.png")), p_correlation, width = 8, height = 5, dpi = 300)


  # 6c. UMAP Visualization (Check if UMAP exists)
  # Determine UMAP reduction name to use
  umap_reduction_name <- config$seurat_input_umap_name %||% "umap"

  if (umap_reduction_name %in% names(seurat_obj@reductions)) {
      message("Generating comparative UMAP plots...")
      top_loading_program <- loading_df %>% arrange(desc(Loading)) %>% slice(1) %>% pull(Program) %>% as.character()
      top_abs_corr_program <- correlation_df %>% arrange(desc(abs(Correlation))) %>% slice(1) %>% pull(Program) %>% as.character()

      p_gene_umap <- FeaturePlot(seurat_obj, features = gene_of_interest, coord.fixed = TRUE, pt.size = 0.1, reduction = umap_reduction_name, order=TRUE) +
                       scale_colour_viridis(option="magma", name=paste(gene_of_interest)) + ggtitle(paste(gene_of_interest, "Expr.")) + theme(plot.title = element_text(hjust = 0.5))
      p_top_loading_umap <- FeaturePlot(seurat_obj, features = top_loading_program, coord.fixed = TRUE, pt.size = 0.1, reduction = umap_reduction_name, order=TRUE) +
                              scale_colour_viridis(option="plasma", name=paste(top_loading_program)) + ggtitle(paste("Highest Loading Prog.")) + theme(plot.title = element_text(hjust = 0.5))
      p_top_corr_umap <- FeaturePlot(seurat_obj, features = top_abs_corr_program, coord.fixed = TRUE, pt.size = 0.1, reduction = umap_reduction_name, order=TRUE) +
                             scale_colour_viridis(option="plasma", name=paste(top_abs_corr_program)) + ggtitle(paste("Highest Corr. Prog.")) + theme(plot.title = element_text(hjust = 0.5))

      combined_umap <- p_gene_umap | p_top_loading_umap | p_top_corr_umap
      print(combined_umap)
      ggsave(file.path(plot_dir, paste0(gene_of_interest, "_nmf_umaps_comparison.png")), combined_umap, width = 15, height = 5, dpi = 300)
  } else {
      message("UMAP reduction '", umap_reduction_name, "' not found, skipping comparative UMAP plots for gene analysis.")
  }


  # 6d. Scatter Plot (Optional)
  message("Generating scatter plot...")
  target_program_for_scatter <- correlation_df %>% arrange(desc(abs(Correlation))) %>% slice(1) %>% pull(Program) %>% as.character()

  if (length(target_program_for_scatter) > 0 && config$celltype_col %in% colnames(seurat_obj@meta.data)) {
      scatter_data <- data.frame(
          GeneExpr = gene_expression,
          ProgramActivity = program_scores[[target_program_for_scatter]],
          Celltype = seurat_obj@meta.data[[config$celltype_col]]
      )

      if (nrow(scatter_data) > config$gene_analysis_scatter_max_cells) {
          set.seed(123)
          scatter_data <- sample_n(scatter_data, config$gene_analysis_scatter_max_cells)
          plot_subtitle <- paste("(Sampled", config$gene_analysis_scatter_max_cells, "cells)")
      } else {
          plot_subtitle <- NULL
      }

      celltypes <- unique(seurat_obj@meta.data[[config$celltype_col]])
      celltype_levels <- levels(factor(seurat_obj@meta.data[[config$celltype_col]]))
      celltype_colors <- scales::hue_pal()(length(celltype_levels))
      names(celltype_colors) <- celltype_levels

      p_scatter <- ggplot(scatter_data, aes(x = ProgramActivity, y = GeneExpr, color = Celltype)) +
          geom_point(alpha = 0.5, size = 0.8) +
          geom_smooth(method = "gam", se = FALSE, color = "black", linetype = "dashed", size = 0.5, aes(group=1)) +
          scale_color_manual(values = celltype_colors, name="Cell Type") +
          theme_bw(base_size = 12) +
          labs(title = paste(gene_of_interest, "Expression vs.", target_program_for_scatter, "Activity"),
               subtitle = plot_subtitle,
               x = paste(target_program_for_scatter, "Activity Score"),
               y = paste(gene_of_interest, "Norm. Expression")) +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
      print(p_scatter)
      ggsave(file.path(plot_dir, paste0(gene_of_interest, "_vs_", target_program_for_scatter, "_scatter.png")), p_scatter, width = 7, height = 6, dpi = 300)
  } else if (!config$celltype_col %in% colnames(seurat_obj@meta.data)) {
      message("Skipping scatter plot as cell type column '", config$celltype_col, "' not found.")
  } else {
       message("No program found with non-zero correlation for scatter plot.")
  }

  message("Specific gene analysis for '", gene_of_interest, "' complete.")
}


# --- 7. Main Workflow Execution ---
run_pipeline <- function(config) {
  message("========================================")
  message("Starting cNMF-hdWGCNA Integration Pipeline")
  message("Time: ", Sys.time())
  message("Output directory: ", normalizePath(config$output_basedir))
  message("========================================")

  load_libraries()
  seurat_obj <- setup_analysis(config)
  saveRDS(seurat_obj, file.path(config$output_basedir, "seurat_object_step1_setup.rds"))

  cnmf_info <- prepare_cnmf_input(seurat_obj, config)

  # --- Determine Optimal K ---
  message("\n>>> PAUSE: Please run the cNMF commands provided above in your terminal (prepare, factorize, combine).")
  readline(prompt="Press [Enter] after external cNMF run is finished...") # Interactive pause

  # Check if user provided K, otherwise attempt auto-selection
  optimal_k <- config$cnmf_optimal_k
  if (is.null(optimal_k)) {
      message("Optimal K not set in config. Attempting automatic selection...")
      optimal_k <- select_optimal_k_auto(cnmf_info$output_dir, cnmf_info$run_name, config$cnmf_k_range)

      if (is.null(optimal_k)) {
          stop("Automatic K selection failed. Please examine the '", config$cnmf_run_name,
               ".k_selection.png' plot in '", normalizePath(cnmf_info$output_dir),
               "' and manually set 'cnmf_optimal_k' in the configuration.")
      } else {
          message("--> Automatically selected K = ", optimal_k)
          message("!!! RECOMMENDATION: Please visually inspect '", config$cnmf_run_name,
                  ".k_selection.png' to verify this choice !!!")
          # Update config for downstream functions if needed (optional)
          # config$cnmf_optimal_k <- optimal_k
      }
  } else {
       message("Using manually specified optimal K = ", optimal_k)
  }
  # --- End Optimal K Determination ---


  # Step 2b: Load cNMF Results using determined optimal_k
  cnmf_results <- load_cnmf_results(seurat_obj, cnmf_info$output_dir, cnmf_info$run_name, optimal_k)
  seurat_obj <- cnmf_results$seurat_obj
  W_matrix <- cnmf_results$W_matrix
  H_matrix <- cnmf_results$H_matrix
  saveRDS(seurat_obj, file.path(config$output_basedir, "seurat_object_step2_cnmf_loaded.rds"))
  saveRDS(W_matrix, file.path(config$output_basedir, "cnmf_W_matrix.rds"))
  saveRDS(H_matrix, file.path(config$output_basedir, "cnmf_H_matrix.rds"))


  # Step 2c: Visualize cNMF Results
  visualize_cnmf_results(seurat_obj, H_matrix, W_matrix, config)

  # Step 3a: Run GSEA on NMF Programs
  gsea_results <- run_gsea_on_nmf(W_matrix, config)
  saveRDS(gsea_results, file.path(config$output_basedir, "gsea_results_list.rds"))

  # Step 3b: Visualize GSEA Results
  visualize_gsea_results(gsea_results, config)

  # Step 4a: Setup hdWGCNA
  seurat_obj <- setup_hdwgcna(seurat_obj, config)
   # --- MANUAL Soft Power SELECTION POINT ---
   if (is.null(config$hdwgcna_soft_power)) {
       message("\n>>> PAUSE: Examine the soft power plots generated in ",
               file.path(config$output_basedir, "plots", "hdwgcna"))
       message(">>> Set 'hdwgcna_soft_power' in the 'config' list below and then continue the script.")
       stop("Please set 'hdwgcna_soft_power' in the configuration before proceeding.")
   }
   message("Continuing with soft power = ", config$hdwgcna_soft_power)
  # --- END MANUAL POINT ---
  saveRDS(seurat_obj, file.path(config$output_basedir, "seurat_object_step4a_hdwgcna_setup.rds"))


  # Step 4b: Run hdWGCNA Network Construction
  seurat_obj <- run_hdwgcna_network(seurat_obj, config)
  saveRDS(seurat_obj, file.path(config$output_basedir, "seurat_object_step4b_hdwgcna_network.rds"))

  # Step 4c: Visualize hdWGCNA Modules
  visualize_hdwgcna_modules(seurat_obj, config)

  # Step 5a: Correlate NMF Programs and WGCNA Modules
  seurat_obj <- correlate_nmf_hdwgcna(seurat_obj, config)
  saveRDS(seurat_obj, file.path(config$output_basedir, "seurat_object_step5a_integration_corr.rds"))
  corr_data <- GetModuleTraitCorrelation(seurat_obj)
  saveRDS(corr_data, file.path(config$output_basedir, "nmf_hdwgcna_correlation_results.rds"))

  # Step 5b: Visualize Module-Trait Correlations
  visualize_module_trait_corr(seurat_obj, config)

  # Step 6: Specific Gene Analysis (if gene specified)
  analyze_gene_nmf_relationship(seurat_obj, W_matrix, config)

  message("========================================")
  message("Pipeline execution finished successfully!")
  message("Time: ", Sys.time())
  message("Final Seurat object saved in: ", config$output_basedir)
  message("========================================")

  return(seurat_obj) # Return the final object
}


# --- Configuration ---
# Define all parameters in a list for easy management
# Load the necessary library for the %||% operator
library(rlang)

config <- list(
  # Input/Output
  input_data_path = "path/to/your/data.rds", # <<< USER: Set path to your input RDS file (SCE or Seurat)
  input_is_seurat_object = FALSE, # <<< USER: Set TRUE if input_data_path points to a pre-processed Seurat object, FALSE if it's an SCE object needing preprocessing.
  output_basedir = "./cnmf_hdwgcna_pipeline_output", # Base directory for all outputs
  celltype_col = "celltype", # <<< USER: Column name in metadata (colData(sce) or seurat_obj@meta.data) with cell type annotations

  # Input Seurat Object Specifics (Only used if input_is_seurat_object = TRUE)
  seurat_input_pca_name = "pca", # <<< USER: Name of the PCA reduction in the input Seurat object (e.g., "pca", "Standardpca")
  seurat_input_umap_name = "umap", # <<< USER: Name of the UMAP reduction in the input Seurat object (e.g., "umap", "StandardpcaUMAP2D")

  # Reticulate/Python (Optional - for saving h5ad or running cNMF via R)
  python_env_path = NULL, # E.g., "my_conda_env_with_cnmf" or "/path/to/venv/bin/python" or NULL
  python_env_type = "conda", # "conda", "virtualenv", or "python_executable"

  # Seurat Preprocessing (Only used if input_is_seurat_object = FALSE)
  seurat_n_hvgs = 2000, # Number of highly variable genes
  seurat_n_pcs = 30, # Number of principal components
  seurat_umap_dims = 30, # Number of PCs to use for UMAP

  # cNMF Parameters
  cnmf_run_name = "scRNA_nmf", # Name for the cNMF run
  cnmf_k_range = seq(5, 20, by = 1), # <<< USER: Set K range to test
  cnmf_n_iter = 100, # Number of NMF iterations per K (recommend >= 100)
  cnmf_seed = 12345,
  cnmf_num_workers = 8, # Number of parallel workers for cNMF factorize step
  cnmf_optimal_k = NULL, # <<< USER: Set this manually based on plot OR leave NULL to attempt automatic selection.

  # cNMF Visualization
  nmf_program_prefix = "NMF_Program_", # Prefix used for NMF program names in metadata/matrices
  cnmf_num_umap_plots = 9, # Max number of NMF program UMAPs to generate
  cnmf_heatmap_ncells = 2000, # Max number of cells for NMF activity heatmap
  cnmf_num_top_genes_heatmap = 15, # Number of top genes per program for W matrix heatmap

  # GSEA Parameters
  gsea_species = "Homo sapiens", # <<< USER: Set species ("Homo sapiens" or "Mus musculus", etc.)
  gsea_categories = "H,C5", # <<< USER: Comma-separated MSigDB categories (e.g., "H", "C5", "C2", "C7")
  gsea_subcategories = "GO:BP,GO:MF", # <<< USER: Comma-separated subcategories if needed (e.g., "GO:BP", "CP:KEGG", "IMMUNESIGDB") or NULL
  gsea_pvalue_cutoff = 0.1, # p-value cutoff for GSEA algorithm
  gsea_min_geneset_size = 10,
  gsea_max_geneset_size = 500,
  # gsea_num_workers = 4, # For parallel execution (requires BiocParallel setup)

  # GSEA Visualization
  gsea_plot_fdr_cutoff = 0.05, # Adjusted p-value cutoff for plotting results
  gsea_plot_top_n_terms = 5, # Max number of top terms per program to show in dot plot

  # hdWGCNA Parameters
  hdwgcna_name = "scNMF_integration", # Name for the hdWGCNA analysis within Seurat object
  hdwgcna_gene_select_method = "variable", # "variable" (uses HVGs), "custom" (provide features), or others
  hdwgcna_custom_features = NULL, # Provide gene vector if method is "custom"
  hdwgcna_assay_for_wgcna = "RNA", # Assay to use for expression data (ensure it exists in the input Seurat object if pre-processed)
  hdwgcna_layer_for_wgcna = "data", # Layer to use ('data' for log-norm, 'scale.data') (ensure it exists)
  hdwgcna_group_name = "all_cells", # Name for the group being analyzed (can be arbitrary if group.by is used)
  hdwgcna_group.by = "celltype", # <<< USER: Metadata column to group cells for network construction (can be NULL)
  hdwgcna_network_type = "signed", # 'signed', 'unsigned', 'signed hybrid'
  hdwgcna_min_rsquared = 0.8, # Minimum R^2 for soft power selection guidance
  hdwgcna_soft_power = NULL, # <<< USER: Set this AFTER running setup_hdwgcna and checking soft power plots
  hdwgcna_tom_name = "NMF_TOM", # Name for the TOM file/matrix
  hdwgcna_compute_module_scores = TRUE, # Whether to compute AddModuleScore for visualization
  hdwgcna_module_score_n_genes = 25, # Number of genes for AddModuleScore

  # hdWGCNA Metacells (Optional)
  hdwgcna_use_metacells = TRUE, # <<< USER: Set TRUE to use metacells, FALSE for single cells
  hdwgcna_metacell_group_by = c("celltype"), # <<< USER: Metadata column(s) to group for metacells
  hdwgcna_metacell_k = 25, # Number of neighbors for metacell construction
  hdwgcna_metacell_max_shared = 12, # Max shared neighbors for metacell construction
  hdwgcna_metacell_ident = "metacell_id", # Name for the new metadata column holding metacell IDs

  # hdWGCNA Visualization
  hdwgcna_num_umap_plots = 9, # Max number of module score UMAPs to generate

  # Integration Visualization
  integration_plot_high_color = "red",
  integration_plot_mid_color = "white",
  integration_plot_low_color = "blue",
  integration_plot_corr_max = 0.8, # Max absolute correlation for color scale

  # Specific Gene Analysis
  analyze_gene_name = NULL, # <<< USER: Set a gene symbol (e.g., "RAPGEF3") or NULL to skip
  gene_analysis_scatter_max_cells = 5000 # Max cells for gene scatter plot
)


# --- Execute Pipeline ---
# It's recommended to run this in two stages due to the manual cNMF and soft power steps.
# Stage 1: Run up to the first pause point.
# Stage 2: Set the required parameters (hdwgcna_soft_power) in the config list above.
#          You can EITHER set cnmf_optimal_k manually OR leave it NULL for auto-selection.
#          Then uncomment and run the main execution function.

# final_seurat_object <- run_pipeline(config)

# Example for interactive execution:
# 1. Set ALL parameters in the `config` list above, including input path and type, and reduction names if needed.
# 2. Run all function definitions (lines before `run_pipeline`).
# 3. Run: `seurat_obj <- setup_analysis(config)`
# 4. Run: `cnmf_info <- prepare_cnmf_input(seurat_obj, config)`
# 5. Execute cNMF commands (prepare, factorize, combine) in terminal.
# 6. Return to R. EITHER set `config$cnmf_optimal_k = your_chosen_k` OR leave it NULL.
# 7. Run the rest of the pipeline using the `run_pipeline` function (uncomment the call)
#    OR run step-by-step starting from the K selection logic inside `run_pipeline`.
#    Remember to set `config$hdwgcna_soft_power` after checking the plot from `setup_hdwgcna`.


config$python_env_path='C:\\Users\\charl\\miniconda3\\envs\\sc\\python.exe'
config$input_data_path='RAPGEF3_cells.rds'
config$input_is_seurat_object=TRUE
config$seurat_input_pca_name='Standardpca'
config$seurat_input_umap_name='StandardUMAP2D'
config$celltype_col='new_celltype'
config$output_basedir='NMF_hdWGCNA'
config$python_env_type='python_executable'
config$hdwgcna_group.by='new_celltype'
run_pipeline(config)