library(Seurat)
library(reticulate)
library(anndata)
library(Biobase)
library(MuSiC)
library(SingleCellExperiment)

########################################
# bulk
data <- read_h5ad("data/pseudo_bulk.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)

expr_matrix <- as.matrix(GetAssayData(data, slot = "counts"))
cell_meta <- data@meta.data

eset <- ExpressionSet(assayData = expr_matrix,
                      phenoData = AnnotatedDataFrame(cell_meta))
label_df <- cell_meta[4:19]

############################################################
# function for deconvolution
deconvolution_music <- function(eset, method, sketch_size){
  # ref
  if (method == "full_data") {
    data <- read_h5ad(paste0("data/", method, ".h5ad" ))
  } else {
    data <- read_h5ad(paste0("experiments/T_and_ILC/", method, "." , sketch_size, ".h5ad" ))
  }
  data$obs_names_make_unique()
  data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)
  sce <- as.SingleCellExperiment(data)

  Est.prop = music_prop(bulk.mtx = exprs(eset), sc.sce = sce, clusters = 'cell_type', #select.ct = select.ct,
                              samples = 'donor_id', verbose = T)
  df2 = Est.prop$Est.prop.weighted
  return(df2)
}

#############################################################
# compute metrics
compute_corr <- function(df1, df2, method) {
  df2 <- as.data.frame(df2)
  all_columns <- colnames(df1)
  missing_cols_df2 <- setdiff(all_columns, colnames(df2))

  for (col in missing_cols_df2) {
    df2[[col]] <- 0
  }

  df1 <- df1[, all_columns, drop = FALSE]
  df2 <- df2[, all_columns, drop = FALSE]
  
  common_rows <- intersect(rownames(df1), rownames(df2))
  df1 <- df1[common_rows, , drop = FALSE]
  df2 <- df2[common_rows, , drop = FALSE]
  
  transposed_df <- as.data.frame(t(df2))
  write.csv(transposed_df, file = paste0("results/", method, ".csv"), row.names = T)

  row_correlations <- sapply(common_rows, function(row_name) {
    vec1 <- as.numeric(df1[row_name, ])
    vec2 <- as.numeric(df2[row_name, ])
    cor(vec1, vec2)
  })
  return(row_correlations)
}

compute_rmse <- function(df1, df2, method) {
  df2 <- as.data.frame(df2)
  all_columns <- colnames(df1)
  missing_cols_df2 <- setdiff(all_columns, colnames(df2))

  for (col in missing_cols_df2) {
    df2[[col]] <- 0
  }

  df1 <- df1[, all_columns, drop = FALSE]
  df2 <- df2[, all_columns, drop = FALSE]
  
  common_rows <- intersect(rownames(df1), rownames(df2))
  df1 <- df1[common_rows, , drop = FALSE]
  df2 <- df2[common_rows, , drop = FALSE]
  

  row_rmses <- sapply(common_rows, function(row_name) {
    vec1 <- as.numeric(df1[row_name, ])
    vec2 <- as.numeric(df2[row_name, ])
    sqrt(mean((vec1 - vec2)^2))
  })
  return(row_rmses)
}

compute_mae <- function(df1, df2, method) {
  df2 <- as.data.frame(df2)
  all_columns <- colnames(df1)
  missing_cols_df2 <- setdiff(all_columns, colnames(df2))

  for (col in missing_cols_df2) {
    df2[[col]] <- 0
  }

  df1 <- df1[, all_columns, drop = FALSE]
  df2 <- df2[, all_columns, drop = FALSE]
  
  common_rows <- intersect(rownames(df1), rownames(df2))
  df1 <- df1[common_rows, , drop = FALSE]
  df2 <- df2[common_rows, , drop = FALSE]

  row_maes <- sapply(common_rows, function(row_name) {
    vec1 <- as.numeric(df1[row_name, ])
    vec2 <- as.numeric(df2[row_name, ])
    mean(abs(vec1 - vec2))
  })
  return(row_maes)
}

eval_method <- function(method, sketch_size) {
  print(method)
  df2 <- deconvolution_music(eset, method, sketch_size)

  corr_list <- compute_corr(label_df, df2, method)
  print('corr: ')
  print(corr_list)
  print(mean(corr_list))
  print(sd(corr_list))
  
  rmse_list <- compute_rmse(label_df, df2, method)
  print('rmse: ')
  print(rmse_list)
  print(mean(rmse_list))
  print(sd(rmse_list))

  mae_list <- compute_mae(label_df, df2, method)
  print('mae: ')
  print(mae_list)
  print(mean(mae_list))
  print(sd(mae_list))
}

###############################################
# deconvolution
eval_method('scValue')
eval_method("Uniform")
eval_method("GeoSketch")
eval_method("Sphetcher")
eval_method("Hopper")
eval_method("KH")
eval_method("scSampler")
eval_method("full_data")
