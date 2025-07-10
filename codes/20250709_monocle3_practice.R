options(future.globals.maxSize = 50.0 * 1e9) ## 50.0 GB
suppressPackageStartupMessages({
  library(rhdf5)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(monocle3)
  library(tibble)
  library(RColorBrewer)
  library(openxlsx)
  library(SeuratWrappers)
  library(pheatmap)
  library(harmony)
  library(sctransform)
  library(glmGamPoi)
  library(future)
  library(reticulate)
  library(future)
  library(future.apply)
  library(igraph)})

f = "/home/hjlee/RCC_upload_final_raw_counts.h5ad"
h5ls(f)
data <- readH5AD(file = f,
                 verbose=T, layers=F, varm=F,
                 obsm=F, varp=F, obsp=F, uns=F)
data = data[, !data$summaryDescription %in% c("Blood", "Fat", "Thrombus")]

# this brings down the number of cells from 270855 to 204255
counts(data) = data@assays@data$X
data@assays@data$X = NULL
data = scater::logNormCounts(data)
data = data[!grepl("^RP", rownames(data)),]
cells_keep <- which(data@colData$percent.mt < 10 & 
                      colSums(data@assays@data@listData[["counts"]]) > 200 & colSums(data@assays@data@listData[["counts"]]) < 5000)
data <- data[, cells_keep]
data <- data[rowSums(data@assays@data@listData[["counts"]]) >= 3, ]
# this brings down the number of cells from 204255 to 182495

data_subset = do.call(cbind, lapply(unique(data$broad_type), function(x) {
subset = data[, data$broad_type == x]
if (ncol(subset) > 15000) {
  set.seed(1)
  subset = subset[, sample(1:ncol(subset), 15000)]
}
return(subset)
}))
data_seurat <- Seurat::as.Seurat(data_subset, counts='counts', data = 'logcounts') # count as raw data; data as log transformed data
data_seurat <- FindVariableFeatures(data_seurat,selection.method = "vst", nfeatures = 2000)
data_seurat <- ScaleData(data_seurat, features = VariableFeatures(object = data_seurat))
data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat)) # select gene targets for Z-normalization
data_seurat <- RunUMAP(data_seurat, dims = 1:10)

DimPlot(data_seurat, reduction = "umap", group.by = "broad_type", label = T, repel = T, raster = T,raster.dpi=c(1024, 1024))
DimPlot(data_seurat, reduction = "umap", group.by = "summaryDescription", label = T, repel = T, raster = T,raster.dpi=c(1024, 1024))
DimPlot(data_seurat, reduction = "umap", group.by = "patient", label = T, repel = T, raster = T,raster.dpi=c(1024, 1024))
DimPlot(data_seurat, reduction = "umap", group.by = "annotation", label = T, repel = T, raster = T,raster.dpi=c(1024, 1024))

# epithelial subset
Epi <- subset(data_seurat, subset = broad_type %in% c("Epi_PT", "Epi_non-PT", "RCC"))
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000) # recalculate HVG
Epi <- ScaleData(Epi, features = rownames(Epi))
Epi <- RunPCA(Epi, features = VariableFeatures(Epi))
Epi <- RunUMAP(Epi, dims = 1:10)

#
DimPlot(Epi, reduction = "umap", group.by = "broad_type", label = T, repel = T)
DimPlot(Epi, reduction = "umap", group.by = "summaryDescription", label = T, repel = T)
DimPlot(Epi, reduction = "umap", group.by = "patient", label = T, repel = T)
DimPlot(Epi, reduction = "umap", group.by = "annotation", label = T, repel = T)

#
data_seurat_harmonized <- readRDS("/home/hjlee/20250630_harmonized_.rds")

DimPlot(data_seurat_harmonized, reduction = "umap", group.by = "broad_type", label = T, raster = T,raster.dpi=c(1024, 1024))

# epithelial subset
Epi_harmonized <- subset(data_seurat_harmonized, 
                         subset = broad_type %in% c("Epi_PT", "Epi_non-PT", "RCC"))
Epi_harmonized <- FindVariableFeatures(Epi_harmonized, selection.method = "vst", nfeatures = 2000) # recalculate HVG
Epi_harmonized <- ScaleData(Epi_harmonized, features = rownames(Epi_harmonized))
Epi_harmonized <- RunPCA(Epi_harmonized, features = VariableFeatures(Epi_harmonized))
Epi_harmonized <- RunUMAP(Epi_harmonized, dims = 1:10, reduction = 'harmony') # treat batch corrected PCA as UMAP input

#
FeaturePlot(Epi_harmonized, features = c("EPCAM", "VIM", "ZEB1"))

DimPlot(Epi_harmonized, reduction = "umap", group.by = "broad_type", label = T, repel = T)
DimPlot(Epi_harmonized, reduction = "umap", group.by = "summaryDescription", label = T, repel = T)
DimPlot(Epi_harmonized, reduction = "umap", group.by = "patient", label = T, repel = T)
DimPlot(Epi_harmonized, reduction = "umap", group.by = "annotation", label = T, repel = T)

DimPlot(Epi_harmonized, reduction = "harmony", group.by = "patient")
DimPlot(Epi_harmonized, reduction = "harmony", group.by = "broad_type")
DimPlot(Epi_harmonized, reduction = "harmony", group.by = "summaryDescription")
DimPlot(Epi_harmonized, reduction = "harmony", group.by = "annotation")

DimPlot(data_seurat_harmonized, reduction = "harmony", group.by = "patient")
DimPlot(data_seurat_harmonized, reduction = "harmony", group.by = "broad_type")
DimPlot(data_seurat_harmonized, reduction = "harmony", group.by = "summaryDescription")
DimPlot(data_seurat_harmonized, reduction = "harmony", group.by = "annotation")

#
cds <- as.cell_data_set(Epi_harmonized, assay = "SCT") # using SeuratWrappers
cds <- estimate_size_factors(cds)
cds_RCC_epi <- cds[, cds@colData@listData[["broad_type"]] %in% c('RCC', 'Epi_PT')]
cds_RCC_epi <- cluster_cells(cds_RCC_epi, reduction_method = "UMAP", resolution = 1e-4)
cds_RCC_epi <- learn_graph(cds_RCC_epi, use_partition = F)
table(cds_RCC_epi@colData@listData[["broad_type"]])

plot_cells(cds_RCC_epi,
           color_cells_by = "broad_type", trajectory_graph_segment_size = 2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           show_trajectory_graph = F,
           label_branch_points=FALSE)
plot_cells(cds_RCC_epi,
           color_cells_by = "broad_type", trajectory_graph_segment_size = 2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)

root_candidates <- which(
  colData(cds_RCC_epi)[, "summaryDescription"] == "Normal kidney" &
    colData(cds_RCC_epi)[, "broad_type"] %in% c("Epi_PT"))
colData(cds_RCC_epi)[['root_or_not']] <- FALSE  
colData(cds_RCC_epi)[['root_or_not']][root_candidates] <- TRUE 

closest_vertex <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex 
closest_vertex <- as.matrix(closest_vertex[colnames(cds_RCC_epi),])
root_pr_nodes <- igraph::V(principal_graph(cds_RCC_epi)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[root_candidates,]))))]
cds_RCC_epi <- order_cells(cds_RCC_epi, root_pr_nodes=root_pr_nodes) # designate starting principal point where pseudo-time begin
cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime <- pseudotime(cds_RCC_epi)
colData(cds_RCC_epi)$pseudotime <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime
colData(cds_RCC_epi)$vertex <- cds_RCC_epi@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]]
table(cds_RCC_epi@colData@listData[["broad_type"]])

root_nodes <- function(cds_RCC_epi, reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                           cds_RCC_epi@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds_RCC_epi@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}
branch_nodes <- function(cds_RCC_epi,reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds_RCC_epi, reduction_method) == FALSE]
  return(branch_points)
}
leaf_nodes <- function(cds_RCC_epi,reduction_method="UMAP"){
  g = principal_graph(cds_RCC_epi)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds_RCC_epi, reduction_method) == FALSE]
  return(leaves)
}
dp_mst <- principal_graph(cds_RCC_epi)[["UMAP"]]
mst_root_nodes <- root_nodes(cds_RCC_epi, "UMAP")
mst_branch_nodes <- branch_nodes(cds_RCC_epi, "UMAP")
mst_leaf_nodes <- leaf_nodes(cds_RCC_epi, 'UMAP')

plot_cells(cds_RCC_epi,
           color_cells_by = "pseudotime", 
           group_cells_by = 'broad_type',
           trajectory_graph_segment_size = 2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = F,
           graph_label_size=1,
           scale_to_range = T)

# cell proportion
cds_sub_df <- data.frame(
  pseudotime = cds_RCC_epi@colData@listData[["pseudotime"]],
  broad_type = cds_RCC_epi@colData@listData[["broad_type"]],
  summaryDescription = cds_RCC_epi@colData@listData[["summaryDescription"]],
  annotation =  cds_RCC_epi@colData@listData[["annotation"]],
  row.names = colnames(cds_RCC_epi))

cds_sub_df_1 <- data.frame(
  pseudotime = cds_RCC_epi@colData@listData[["pseudotime"]],
  broad_type = cds_RCC_epi@colData@listData[["broad_type"]],
  summaryDescription = cds_RCC_epi@colData@listData[["summaryDescription"]],
  annotation =  cds_RCC_epi@colData@listData[["annotation"]],
  vertex = cds_RCC_epi@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]],
  row.names = colnames(cds_RCC_epi))

closest_vertex <- cds_RCC_epi@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]]
names(closest_vertex) <- colnames(cds_RCC_epi)

cds_sub_df$pseudo_bin <- cut(cds_RCC_epi$pseudotime, breaks = 8)
cds_sub_df_1$pseudo_bin <- cut(cds_RCC_epi$pseudotime, breaks = 8)

plot_df <- cds_sub_df_1 %>%
  group_by(pseudo_bin, annotation) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(pseudo_bin) %>%
  mutate(freq = n / sum(n))
plot_df$annotation <- as.factor(plot_df$annotation)

n_colors <- length(unique(plot_df$annotation))  
my_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n_colors)

ggplot(plot_df, aes(x = pseudo_bin, y = freq, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Pseudotime bin", y = "Proportion", fill = "Annotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# terminal node test
getwd()
for(i in 1:5){
  x <- mst_leaf_nodes[i]
  ENDnodes <- c(x)
  path.nodes <- NULL
  for (end in ENDnodes){
    path <- igraph::shortest_paths(
      dp_mst, from <- root_pr_nodes, to = end, mode = "all",
      algorithm = "unweighted")
    nodes <- path$vpath[[1]]
    path.nodes <- c(path.nodes, nodes)
  }
  path.nodes <- unique(path.nodes)
  cells.branch <- closest_vertex[closest_vertex %in% path.nodes,]
  cds_RCC_epi@colData$sub_trajectory <- names(cds_RCC_epi@clusters$UMAP$clusters) %in% names(cells.branch)
  p1 <- plot_cells(cds_RCC_epi, 
                   color_cells_by = "sub_trajectory", trajectory_graph_segment_size = 1,
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=3) +
    ggtitle(paste0('mst_leaf_nodes_',x,'_result')) +
    theme(plot.title = element_text(hjust = 0.1, size=12, face = 'bold'))
  jpeg(filename = paste0('trajectory',x,'_leafnode_results','.jpeg'), width = 10, height = 7, units = 'cm', res = 300)
  plot(p1)
  dev.off()
}
# designate coordinate for branch
ENDnodes <- c("Y_88")
path.nodes <- NULL
for (end in ENDnodes){
  path <- igraph::shortest_paths(
    dp_mst, from <- root_pr_nodes, to = end, mode = "all",
    algorithm = "unweighted")
  nodes <- path$vpath[[1]]
  path.nodes <- c(path.nodes, nodes)
}
path.nodes <- unique(path.nodes)
cells.branch <- closest_vertex[closest_vertex %in% path.nodes,] # select cell that included in path.node
cds_RCC_epi@colData$sub_trajectory <- names(cds_RCC_epi@clusters$UMAP$clusters) %in% names(cells.branch) # whether the cell 
cds_RCC_epi@colData$optimal_path <- rownames(cds_RCC_epi@colData) %in% names(cells.branch)

# differential expression gene across a single-cell trajectory 
cds_sub <- cds_RCC_epi[, cds_RCC_epi@colData$sub_trajectory]
table(cds_sub@colData@listData[["broad_type"]])

deg_res <- graph_test(cds_sub, neighbor_graph="principal_graph", cores = 4) # identify genes with interesting patterns of expression that fall only within the region of the trajectory
deg_res_sig <- rownames(subset(deg_res, morans_I > 0 & q_value < 0.05))
pr_deg_ids_1 <- rownames(subset(deg_res, morans_I > 0))

cds_sub <- preprocess_cds(cds_sub, num_dim = 10) # PCA re 
cds_sub <- reduce_dimension(cds_sub) # UMAP re 

gene_module_df <- find_gene_modules(cds_sub, resolution=0.1, random_seed = 123) # 
gene_module_df_sig <- gene_module_df %>% filter(id %in% deg_res_sig)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_sub)), 
                                cell_group=cds_sub@colData@listData[["annotation"]])

agg_mat <- aggregate_gene_expression(cds_sub, gene_module_df_sig, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

my_palette <- colorRampPalette(c("steelblue3", "white", "tomato3"))(1200)
min_val <- min(agg_mat)
max_val <- max(agg_mat)
mid_val <- 0
breaks <- c(
  seq(min_val, mid_val, length.out = 600),
  seq(mid_val, max_val, length.out = 601)[-1])
pheatmap::pheatmap(t(agg_mat), angle_col = 90, cutree_cols = 3,
                   , cutree_rows = 2, color = my_palette, breaks = breaks, 
                   scale="row", clustering_method="ward.D2", labels_col = F, width = 7, height = 10)