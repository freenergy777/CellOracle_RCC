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
})

data_seurat <- readRDS("/home/hjlee/20250420_3k_SCTransform_suerat.rds")
data_seurat_umap_coords <- as.data.frame(data_seurat@reductions[["umap"]]@cell.embeddings) 

keep_cells <- rownames(data_seurat_umap_coords)[
  data_seurat_umap_coords[,1] >= -10 &
    data_seurat_umap_coords[,1] <= 3.5 &
    data_seurat_umap_coords[,2] >= 0]

data_seurat <- subset(data_seurat, cells = keep_cells)
table(data_seurat@meta.data[["broad_type"]])

cds <- as.cell_data_set(data_seurat, assay = "SCT") # using SeuratWrappers
cds <- estimate_size_factors(cds)
cds_RCC_epi <- cds[, cds@colData@listData[["broad_type"]] %in% c('RCC', 'Epi_PT', 'Epi_non-PT')]
cds_RCC_epi <- cluster_cells(cds_RCC_epi, reduction_method = "UMAP", resolution = 1e-3)
cds_RCC_epi <- learn_graph(cds_RCC_epi)
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
    colData(cds_RCC_epi)[, "broad_type"] %in% c("Epi_PT")
)

colData(cds_RCC_epi)[['root_or_not']] <- FALSE  
colData(cds_RCC_epi)[['root_or_not']][root_candidates] <- TRUE 

closest_vertex <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex 
closest_vertex <- as.matrix(closest_vertex[colnames(cds_RCC_epi),])
root_pr_nodes <- igraph::V(principal_graph(cds_RCC_epi)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[root_candidates,]))))]
cds_RCC_epi <- order_cells(cds_RCC_epi, root_pr_nodes=root_pr_nodes) # designate starting principal point where pseudo-time begin
cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime <- pseudotime(cds_RCC_epi)
colData(cds_RCC_epi)$pseudotime <- cds_RCC_epi@principal_graph_aux[["UMAP"]]$pseudotime
table(cds_RCC_epi@colData@listData[["broad_type"]])

#
wanted_gene <- 'VIM'
gene_index <- which(rownames(cds_RCC_epi) == wanted_gene)
gene_counts_VIM <- cds_RCC_epi@assays@data@listData[["counts"]][gene_index, ]
summary(gene_counts_VIM)
gene_logcounts_VIM <- cds_RCC_epi@assays@data@listData[["logcounts"]][gene_index, ]
summary(gene_logcounts_VIM)

#
features1 <- c("NDUFA4L2","CA9","KRT18","KRT8","KLRD1","KLRB1","GNLY","NKG7","CD3D","CD3E","CD8A","CD8B","IL7R","CD68","GPNMB","SLC40A1","ACTA2","PDGFRB","PECAM1","KDR","CDH5","LYZ","KCNQ1OT1","VWF","KRT19","WFDC2")
DotPlot(data_seurat, features = features1, cols = c("steelblue3", "tomato3"), dot.scale = 8, dot.min = 0.1) + RotatedAxis()

# choose sub trajectory manually
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
           color_cells_by = "pseudotime", trajectory_graph_segment_size = 2,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = F,
           graph_label_size=3)

for(i in 1:20){
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
  jpeg(filename = paste0('trajectory',x,'_leafnode_results','.jpeg'), width = 12, height = 7, units = 'cm', res = 300)
  plot(p1)
  dev.off()
}

# designate Y_13 for branch
ENDnodes <- c("Y_402")
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
write.xlsx(deg_res, 'D:/2025.05/monocle3_practice_2/20250522_Y_249_deg_res_sub_trajectory.xlsx', rowNames=T)
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

cds_sub_df <- data.frame(
  pseudotime = cds_sub@colData@listData[["pseudotime"]],
  broad_type = cds_sub@colData@listData[["broad_type"]],
  summaryDescription = cds_sub@colData@listData[["summaryDescription"]],
  annotation =  cds_sub@colData@listData[["annotation"]],
  row.names = colnames(cds_sub))

cds_sub_df_1 <- data.frame(
  pseudotime = cds_sub@colData@listData[["pseudotime"]],
  cell_type = cds_sub@colData@listData[["broad_type"]],
  annotation = cds_sub@colData@listData[["annotation"]], 
  row.names = colnames(cds_sub))

cds_sub_df$pseudo_bin <- cut(cds_sub_df$pseudotime, breaks = 10)

plot_df <- cds_sub_df %>%
  group_by(pseudo_bin, annotation) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(pseudo_bin) %>%
  mutate(freq = n / sum(n))

n_colors <- length(unique(plot_df$annotation))  
my_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n_colors)

ggplot(plot_df, aes(x = pseudo_bin, y = freq, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Pseudotime bin", y = "Proportion", fill = "Annotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Moran I index
genes_to_plot <- rownames(deg_res)[order(deg_res$morans_I, decreasing = TRUE)][1:50]
exprs_mat <- logcounts(cds_sub)[genes_to_plot, ]
summary(exprs_mat@x)
clustering <- kmeans(exprs_mat, centers = 5)
exprs_ordered <- exprs_mat[, order(cds_sub@colData$pseudotime, na.last = NA)] # make pseudotime order
gene_order <- names(sort(clustering$cluster))  # gene names
exprs_final <- exprs_ordered[gene_order, ]

label_col <- rep("", ncol(exprs_final))
label_col[which(colnames(exprs_final) == "5739STDY8351218_ACGTCAAGTTGGGACA-1")] <- "RCC"
label_col[which(colnames(exprs_final) == "5739STDY8351266_GTAGTCAGTGTAACGG-1")] <- "Epi_non-PT"

pheatmap(exprs_final,
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         labels_col = label_col,
         annotation_row = data.frame(Cluster = factor(clustering$cluster[gene_order])))

# 
marker_genes <- c('NDUFA4L2', 'VIM', 'DEFB1', 'WFDC2', 'APOE', 'GATM')
marker_genes_RCC <- c('NDUFA4L2', 'VIM')
marker_genes_epi_non <- c('DEFB1', 'WFDC2')
marker_genes_epi <- c('APOE', 'GATM')

AFD_lineage_cds <- cds_RCC_epi[cds_RCC_epi@rowRanges@partitioning@NAMES %in% marker_genes,
                               cds_RCC_epi@colData@listData[["sub_trajectory"]]]
AFD_lineage_cds <- cluster_cells(AFD_lineage_cds, reduction_method = "UMAP")
AFD_lineage_cds <- learn_graph(AFD_lineage_cds, use_partition = T)
AFD_lineage_cds <- order_cells(AFD_lineage_cds)

AFD_lineage_cds$pseudo_bin <- cut(AFD_lineage_cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]], breaks = 10)

rowData(AFD_lineage_cds)$gene_short_name <- rownames(AFD_lineage_cds)

plot_genes_in_pseudotime(AFD_lineage_cds, 
                         min_expr = 1e-5,
                         color_cells_by="pseudo_bin",
                         cell_size = 0.1) +
  scale_y_continuous(name = "Normalized expression",
                     limits = c(0,75)) 
# gene order correction
AFD_marker_cds <- AFD_lineage_cds[rowData(AFD_lineage_cds)$gene_short_name %in% marker_genes, ]
rowData(AFD_marker_cds)$gene_short_name <- factor(rowData(AFD_marker_cds)$gene_short_name, levels = marker_genes)

plot_genes_in_pseudotime(AFD_marker_cds,
                         panel_order = marker_genes,
                         min_expr = 1e-5,
                         color_cells_by="pseudo_bin",
                         cell_size = 0.1) +
  scale_y_continuous(name = "Normalized expression",
                     limits = c(0,75)) 

# extract pseudotime and broad_type label
cds_sub_df <- data.frame(
  pseudotime = cds_sub@colData@listData[["pseudotime"]],
  broad_type = cds_sub@colData@listData[["broad_type"]],
  row.names = colnames(cds_sub))
cds_sub_df <- cds_sub_df[order(cds_sub_df$pseudotime),]
write.csv(cds_sub_df, "/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250625_order.csv")

#
f = "/home/hjlee/RCC_upload_final_raw_counts.h5ad"
h5ls(f)

data <- readH5AD(file = f,
                 verbose=T, layers=F, varm=F,
                 obsm=F, varp=F, obsp=F, uns=F)
data = data[, !data$summaryDescription %in% c("Blood", "Fat", "Thrombus")] # this brings down the number of cells from 270855 to 204255

counts(data) = data@assays@data$X
data@assays@data$X = NULL
data = data[!grepl("^RP", rownames(data)),]
data <- data[,data@colData@listData[["percent.mt"]]<=10]
data$tissue_major = paste0(data$summaryDescription, "__", 
                           data$broad_type)
keep = names(which(table(paste0(data$summaryDescription, "__", data$broad_type)) > 100))
data = data[, data$tissue_major %in% keep]
data <- data[rownames(data) %in% deg_res_sig,
             colnames(data) %in% rownames(cds_sub_df)]
data_rawcount <- as.data.frame(data@assays@data@listData[["counts"]])
write.csv(data_rawcount, file = '/home/hjlee/monocle3,pyscenic,palantir,iqcell/20250625_renal_raw.csv')
