library(Signac)
library(Seurat)
library(future.apply)
library(stringr)

plan(multisession, workers = 8)

input_dir <- "/media/AnalysisDisk1/xuzekai/20250715_BMI/ATAC/seurat/"
output_dir <- "/media/AnalysisDisk1/xuzekai/20250715_BMI/ATAC/peak_change_LR_TFIDF/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

cell_types_all <- str_replace(basename(rds_files), "\\.rds$", "")

csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
existing_cell_types <- str_replace(csv_files, "\\.csv$", "")

files_to_process <- rds_files[basename(rds_files) %in% paste0(existing_cell_types, ".rds")]

min_pct <- 0.01
logfc_threshold <- 0.1

if (length(files_to_process) > 0) {
  message("开始处理 ", length(files_to_process), " 个新细胞类型...")
  results <- future_lapply(files_to_process, function(file_path) {
    
    cell_type <- str_replace(basename(file_path), "\\.rds$", "")
    output_file <- file.path(output_dir, paste0(cell_type, ".csv"))
    
    message(" Processing: ", cell_type)
    
    tryCatch({
      seurat <- readRDS(file_path)
      
      if (!"BMI_group" %in% colnames(seurat@meta.data)) {
        warning(cell_type, " missing BMI_group, skipped.")
        return(NULL)
      }
      
      bmi_groups <- unique(seurat$BMI_group)
      bmi_groups <- bmi_groups[!is.na(bmi_groups)]
      if (length(bmi_groups) < 2) {
        warning(cell_type, " has <2 BMI groups, skipped.")
        return(NULL)
      }
      
      group_pairs <- combn(bmi_groups, 2, simplify = FALSE)
      all_markers <- list()
      
      for (pair in group_pairs) {
        ident_1 <- pair[1]
        ident_2 <- pair[2]
        comparison <- paste(ident_1, ident_2, sep = "_vs_")
        
        message("  Comparing: ", comparison)
        
        seurat_subset <- subset(seurat, subset = BMI_group %in% c(ident_1, ident_2))
        Idents(seurat_subset) <- seurat_subset$BMI_group
        seurat_subset <- RunTFIDF(seurat_subset)
        
        markers <- FindMarkers(
          seurat_subset,
          slot = "data",
          ident.1 = ident_1,
          ident.2 = ident_2,
          test.use = "LR",
          min.pct = min_pct,
          logfc_threshold = logfc_threshold,
          latent.vars = "nCount_RNA",
          random.seed = 1
        )
        
        markers$comparison <- comparison
        markers$gene <- rownames(markers)
        markers$cell_type <- cell_type
        
        all_markers[[comparison]] <- markers
      }
      
      if (length(all_markers) > 0) {
        combined_markers <- do.call(rbind, all_markers)
        rownames(combined_markers) <- NULL
        write.csv(combined_markers, output_file, row.names = FALSE)
        message("[", Sys.time(), "] Saved: ", output_file)
      }
      
      return(output_file)
      
    }, error = function(e) {
      warning("Error processing ", cell_type, ": ", e$message)
      return(NULL)
    })
  })
} else {
  message("所有细胞类型已处理完毕")
}

print("All cell types processed with multiprocessing.")