library(Seurat)
library(ggplot2)
library(patchwork)

colors <- c('scValue-core' = '#FF4500',
            'scValue' = '#FB61D7', 
            'Uniform' = '#8FB4DC', 
            'GeoSketch' = '#FFDD8E',  
            'Sphetcher' = '#7AC3DF', 
            'Hopper' = '#AC99D2', 
            'KH' = '#EB7E60', 
            'scSampler' = '#70CDBE',
            'Full data' = '#BEBEBE')

##########################
# find top3 markers for each cell type
dat <- readRDS('rds_data/full_data.rds')
dat[["RNA"]]$data <- dat[["RNA"]]$counts

cell_types <- unique(dat$CellType)
Idents(dat) <- dat$CellType

marker_top3_list <- list()

for (cell_type in cell_types) {
  
  if (cell_type == 'Unassigned') next
  
  marker_df <- read.csv(paste0('cell_type_specific_markers/', cell_type, '.csv'))
  marker_list <- marker_df$marker
  marker_exist_list <- intersect(marker_list, rownames(dat))
  
  markers_de <- FindMarkers(
    object = dat,
    ident.1 = cell_type,
    ident.2 = NULL,
    features = marker_exist_list,
    only.pos = TRUE
  )
  
  markers_top3 <- row.names(markers_de)[1:3]
  marker_top3_list[[cell_type]] <- markers_top3
  print(cell_type)
  print(markers_top3)
}

################################
# merge datasets
dat <- readRDS('rds_data/full_data.rds')
dat$dataset <- 'Full data'
dat_scValue_core <- readRDS('rds_data/scValue-core.rds')
dat_scValue_core$dataset <- 'scValue-core'
dat_scValue <- readRDS('rds_data/scValue.rds')
dat_scValue$dataset <- 'scValue'
dat_Uniform <- readRDS('rds_data/Uniform.rds')
dat_Uniform$dataset <- 'Uniform'
dat_GeoSketch <- readRDS('rds_data/GeoSketch.rds')
dat_GeoSketch$dataset <- 'GeoSketch'
dat_Sphetcher <- readRDS('rds_data/Sphetcher.rds')
dat_Sphetcher$dataset <- 'Sphetcher'
dat_Hopper <- readRDS('rds_data/Hopper.rds')
dat_Hopper$dataset <- 'Hopper'
dat_KH <- readRDS('rds_data/KH.rds')
dat_KH$dataset <- 'KH'
dat_scSampler <- readRDS('rds_data/scSampler.rds')
dat_scSampler$dataset <- 'scSampler'

datasets <- c("Full data", "scValue-core", "scValue",
              "Uniform", "GeoSketch", "Sphetcher",
              "Hopper", "KH", "scSampler")

dat_merged <- merge(dat, y = c(dat_scValue_core,
                               dat_scValue,
                               dat_Uniform,
                               dat_GeoSketch, 
                               dat_Sphetcher,
                               dat_Hopper,
                               dat_KH,
                               dat_scSampler),
                    add.cell.ids = datasets)

Idents(dat_merged) <- dat_merged$CellType

################################################
# violin plot of top 3 markers for each cell type

for (cell_type in names(marker_top3_list)) {
  
  print(paste('Plotting', cell_type))
  
  dat_cell_type <- subset(dat_merged, idents = cell_type)
  dat_cell_type$dataset <- factor(dat_cell_type$dataset, levels = datasets)
  
  markers_top3 <- marker_top3_list[[cell_type]]
  
  p1 <- VlnPlot(dat_cell_type, features = markers_top3[1], group.by = "dataset", pt.size = 0, cols = colors) +
    ggtitle(markers_top3[1])+
    theme(axis.title.x = element_blank())
  
  p2 <- VlnPlot(dat_cell_type, features = markers_top3[2], group.by = "dataset", pt.size = 0, cols = colors) +
    ggtitle(markers_top3[2])+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p3 <- VlnPlot(dat_cell_type, features = markers_top3[2], group.by = "dataset", pt.size = 0, cols = colors) +
    ggtitle(markers_top3[3])+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(
      title = cell_type,
      theme = theme(
        plot.title = element_text(size = 14, 
                                  face = "bold",
                                  hjust = 0, 
                                  margin = margin(b = 10))
      )
    )
  #print(combined_plot)
  
  ggsave(
    filename = paste0(cell_type, '.tiff'),   
    plot     = combined_plot,         
    width    = 12,                
    height   = 6,                 
    dpi      = 300,     
    bg = 'white'
  )
}


