library(SpatialExperiment)
library(SpatialArtifacts)
library(arrow)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(scales)
library(terra)
library(DropletUtils)

zip_file_path <- "/users/hjialihe/.cache/R/BiocFileCache/10049d62745472_VisiumHD_HumanColon_Oliveira.zip"

unzip_dir <- tempfile()
unzip(zip_file_path, exdir = unzip_dir)

path_to_16um_data <- file.path(unzip_dir, "binned_outputs", "square_016um")
h5_file_16um      <- file.path(path_to_16um_data, "filtered_feature_bc_matrix.h5")
spatial_dir_16um  <- file.path(path_to_16um_data, "spatial")
spatial_file_16um <- file.path(spatial_dir_16um, "tissue_positions.parquet")

h5_files <- list.files(unzip_dir, pattern = "filtered_feature_bc_matrix.h5", recursive = TRUE, full.names = TRUE)
h5_target <- h5_files[grep("square_016um", h5_files)]
if (length(h5_target) == 0) h5_target <- h5_files[1]

parquet_files <- list.files(unzip_dir, pattern = "tissue_positions.parquet", recursive = TRUE, full.names = TRUE)
parquet_target <- parquet_files[grep("square_016um", parquet_files)]
if (length(parquet_target) == 0) parquet_target <- parquet_files[1]

sce_temp <- read10xCounts(h5_target, col.names = TRUE)

xyz <- arrow::read_parquet(parquet_target)
xyz_df <- as.data.frame(xyz)
rownames(xyz_df) <- xyz_df$barcode
xyz_matched <- xyz_df[colnames(sce_temp), ]

coords_matrix <- as.matrix(xyz_matched[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
rownames(coords_matrix) <- colnames(sce_temp)

spe_hd <- SpatialExperiment(
  assays = assays(sce_temp),
  rowData = rowData(sce_temp),
  colData = DataFrame(xyz_matched),
  spatialCoords = coords_matrix
)

spe <- spe_hd
spe$in_tissue <- as.logical(spe$in_tissue)
spe$sample_id <- "VisiumHD_HumanColon_Oliveira"
colData(spe)$pxl_col_in_fullres <- spatialCoords(spe)[, 1]
colData(spe)$pxl_row_in_fullres <- spatialCoords(spe)[, 2]

spe <- scuttle::addPerCellQCMetrics(spe)

if ("sum" %in% colnames(colData(spe))) {
  colnames(colData(spe))[colnames(colData(spe)) == "sum"] <- "sum_umi"
}

spe_detected <- detectEdgeArtifacts(
  spe                  = spe,
  platform             = "visiumhd",
  resolution           = "16um",
  qc_metric            = "sum_umi",
  samples              = "sample_id",
  mad_threshold        = 3,
  col_x                = "array_col",
  col_y                = "array_row",
  buffer_width_um      = 160,
  min_cluster_area_um2 = 1280,
  name                 = "edge_dryspot_hd",
  keep_intermediate    = FALSE,
  verbose              = TRUE
)

spe_classified <- classifyEdgeArtifacts(
  spe       = spe_detected,
  qc_metric = "sum_umi",
  samples   = "sample_id",
  min_spots = 200,
  name      = "edge_dryspot_hd"
)

spe_classified <- spe_classified[, colSums(counts(spe_classified)) > 0]

set.seed(123)
spe_classified <- logNormCounts(spe_classified)
spe_classified <- runPCA(spe_classified, ncomponents = 30)
spe_classified <- runUMAP(spe_classified)

g <- buildSNNGraph(spe_classified, k = 20, use.dimred = 'PCA')
clust <- igraph::cluster_louvain(g)$membership
spe_classified$cluster <- factor(clust)

df_col <- as.data.frame(colData(spe_classified))
df_umap <- as.data.frame(reducedDim(spe_classified, "UMAP"))
colnames(df_umap) <- c("UMAP1", "UMAP2")
plot_data_all <- cbind(df_col, df_umap)

if (!"pxl_col_in_fullres" %in% colnames(plot_data_all)) {
  plot_data_all <- cbind(plot_data_all, as.data.frame(spatialCoords(spe_classified)))
}

plot_data_all <- plot_data_all[plot_data_all$in_tissue, ]

cls_col <- "edge_dryspot_hd_classification"
plot_data_all$class_factor <- factor(plot_data_all[[cls_col]], levels = c(
  "large_edge_artifact", "small_edge_artifact",
  "large_interior_artifact", "small_interior_artifact", "not_artifact"
))
levels(plot_data_all$class_factor) <- c(
  "Large Edge", "Small Edge", "Large Interior", "Small Interior", "Normal Spot"
)

artifact_colors <- c(
  "Normal Spot"    = "lightgray",
  "Small Edge"     = "#FFB347",
  "Large Edge"     = "#FF4500",
  "Small Interior" = "#00CED1",
  "Large Interior" = "#0047AB"
)

n_clusters <- length(unique(plot_data_all$cluster))
cluster_colors <- scales::hue_pal()(n_clusters)
names(cluster_colors) <- levels(plot_data_all$cluster)

std_threshold <- 500
plot_data_all$std_qc_status <- ifelse(plot_data_all$sum_umi < std_threshold, "Fail", "Pass")

hd_point_size <- 0.8
hd_alpha <- 1
hd_stroke <- 0
x_limits <- range(plot_data_all$pxl_col_in_fullres)
y_limits <- range(plot_data_all$pxl_row_in_fullres)

my_theme <- theme_void() +
  theme(plot.title = element_text(size=16, face="bold", hjust=0),
        legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11),
        plot.margin = margin(2, 2, 2, 2))

theme_vio <- theme_classic() +
  theme(plot.title = element_text(size=16, face="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.title.x = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

p_a <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=log10(sum_umi))) +
  geom_point(size=hd_point_size, stroke=0, alpha=hd_alpha) +
  scale_color_viridis_c(option="magma", name="log10(UMI)") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle("A. Total UMI Distribution") + my_theme + theme(legend.position="right")

p_b <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=std_qc_status)) +
  geom_point(size=hd_point_size, stroke=0, alpha=hd_alpha) +
  scale_color_manual(values=c("Fail"="black", "Pass"="lightgray"), name="Std QC") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle(paste0("B. Standard QC (UMI < ", std_threshold, ")")) + my_theme

plot_data_sorted <- plot_data_all[order(ifelse(plot_data_all$class_factor == "Normal Spot", 1, 2)), ]
p_c <- ggplot(plot_data_sorted, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=class_factor)) +
  geom_point(size=hd_point_size, stroke=0, alpha=hd_alpha) +
  scale_color_manual(values=artifact_colors, name="Type") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle("C. SpatialArtifacts Classification") + my_theme + theme(legend.position="right")

target_cluster_id <- "1"
cluster1_data <- subset(plot_data_all, cluster == target_cluster_id)
p_d <- ggplot() +
  geom_point(data=plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres),
             color="grey95", size=hd_point_size, stroke=0) +
  geom_point(data=cluster1_data, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=class_factor),
             size=hd_point_size, stroke=0, alpha=1) +
  scale_color_manual(values=artifact_colors, name="Type") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle(paste0("D. Anatomy of Cluster ", target_cluster_id)) +
  my_theme + theme(legend.position="right")

p_e <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors) +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  ggtitle("E. Artifact Composition") + ylab("Percentage") + xlab("Cluster ID") +
  theme(plot.title=element_text(size=16, face="bold"),
        legend.position="none",
        axis.text.x=element_text(size=10, angle=90),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        plot.margin=margin(2, 2, 2, 2))

global_cutoff <- 500
plot_data_all$is_global_low <- plot_data_all$sum_umi < global_cutoff
plot_data_global <- plot_data_all[order(plot_data_all$is_global_low), ]

p_f <- ggplot(plot_data_global, aes(x=UMAP1, y=UMAP2, color=is_global_low)) +
  geom_point(size=0.5, alpha=0.6, stroke=0) +
  scale_color_manual(values=c("FALSE"="#E0E0E0", "TRUE"="#D32F2F"),
                     labels=c("Pass", paste0("Fail (<", global_cutoff, ")"))) +
  theme_classic() +
  ggtitle(paste0("F. Global Threshold (UMI < ", global_cutoff, ")")) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        plot.title=element_text(size=16, face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        plot.margin=margin(2, 2, 2, 2))

p_g <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  scale_fill_manual(values=artifact_colors) +
  ggtitle("G. UMI") + ylab("log10(UMI)") + theme_vio

p_h <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  scale_fill_manual(values=artifact_colors) +
  ggtitle("H. Genes") + ylab("log10(Genes)") + theme_vio

final_plot_complete <- (p_a | p_b | p_c) /
                       (p_d | p_e | p_f) /
                       (p_g | p_h) +
                       plot_layout(heights=c(1.2, 1.0, 0.9)) &
                       theme(plot.margin=margin(1, 1, 1, 1))

ggsave("Figure4_VisiumHD_Colon.png", plot=final_plot_complete, width=15, height=16, dpi=300, bg="white")

# Supplementary Figure 4: Comprehensive Validation
n_clust <- length(unique(plot_data_all$cluster))
cols_bio <- scales::hue_pal()(n_clust)
names(cols_bio) <- levels(plot_data_all$cluster)
 
cols_hl <- rep("#E0E0E0", n_clust)
names(cols_hl) <- levels(plot_data_all$cluster)
cols_hl[target_cluster_id] <- "#FF4500"
 
hd_shape <- 15
hd_size_square <- 0.8
 
my_theme_supp <- theme_void() +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5))
 
df_clean <- plot_data_all[plot_data_all$class_factor == "Normal Spot", ]
 
p_s1 <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(shape=hd_shape, size=hd_size_square, stroke=0) +
  scale_color_manual(values=cols_bio) +
  coord_fixed() + scale_y_reverse() + my_theme_supp +
  ggtitle("S1. Colorectal (OSTA) Visium HD\nSpatial Clusters (All)") +
  theme(legend.position="none")
 
p_s2 <- ggplot(df_clean, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(shape=hd_shape, size=hd_size_square, stroke=0) +
  scale_color_manual(values=cols_bio) +
  coord_fixed() + scale_y_reverse() + my_theme_supp +
  ggtitle("S2. Colorectal (OSTA) Visium HD\nSpatial Clusters (Filtered)") +
  theme(legend.position="none")
 
df_ord <- plot_data_all[order(ifelse(plot_data_all$cluster == target_cluster_id, 2, 1)), ]
 
p_u1 <- ggplot(df_ord, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, stroke=0, alpha=0.8) +
  scale_color_manual(values=cols_hl) +
  theme_classic() +
  ggtitle(paste0("U1. Colorectal (OSTA) Visium HD\nUMAP Before (Highlight C", target_cluster_id, ")")) +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5),
        axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
 
p_u2 <- ggplot(df_clean, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, stroke=0, alpha=0.6) +
  scale_color_manual(values=cols_hl) +
  theme_classic() +
  ggtitle("U2. Colorectal (OSTA) Visium HD\nUMAP After (Filtered)") +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5),
        axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
 
run_pb_pca <- function(spe_obj, cluster_ids) {
  pb <- aggregateAcrossCells(spe_obj, ids=cluster_ids)
  pb <- pb[, colSums(counts(pb)) > 0, drop=FALSE]
  sizeFactors(pb) <- colSums(counts(pb))
  pb <- logNormCounts(pb)
  pb <- runPCA(pb, ncomponents=min(2, ncol(pb)-1))
  df <- as.data.frame(reducedDim(pb, "PCA"))
  df$cluster <- colnames(pb)
  var_exp <- attr(reducedDim(pb, "PCA"), "percentVar")
  return(list(df=df, pc1_var=round(var_exp[1], 1), pc2_var=round(var_exp[2], 1)))
}
 
pb_res_before <- run_pb_pca(spe_classified, spe_classified$cluster)
pb_df_before <- pb_res_before$df
counts_all <- table(spe_classified$cluster)
pct_all <- round((counts_all / sum(counts_all)) * 100, 1)
pb_df_before$cluster_label <- paste0(pb_df_before$cluster, "\n(", pct_all[as.character(pb_df_before$cluster)], "%)")
 
p_p1 <- ggplot(pb_df_before, aes(x=PC1, y=PC2, color=cluster, label=cluster_label)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(color="black", size=4, fontface="bold", box.padding=0.6,
                  point.padding=0.5, min.segment.length=0, max.overlaps=Inf) +
  scale_color_manual(values=cols_hl) + theme_classic() +
  ggtitle("P1. Colorectal (OSTA) Visium HD\nPseudo-bulk PCA (Before)") +
  xlab(paste0("PC1 (", pb_res_before$pc1_var, "%)")) +
  ylab(paste0("PC2 (", pb_res_before$pc2_var, "%)")) +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5)) +
  scale_x_continuous(expand=expansion(mult=0.2)) +
  scale_y_continuous(expand=expansion(mult=0.2))
 
spe_clean_obj <- spe_classified[, spe_classified$edge_dryspot_hd_classification == "not_artifact"]
spe_clean_obj$cluster <- droplevels(spe_clean_obj$cluster)
 
pb_res_after <- run_pb_pca(spe_clean_obj, spe_clean_obj$cluster)
pb_df_after <- pb_res_after$df
counts_clean <- table(spe_clean_obj$cluster)
pct_clean <- round((counts_clean / sum(counts_clean)) * 100, 1)
pb_df_after$cluster_label <- paste0(pb_df_after$cluster, "\n(", pct_clean[as.character(pb_df_after$cluster)], "%)")
 
p_p2 <- ggplot(pb_df_after, aes(x=PC1, y=PC2, color=cluster, label=cluster_label)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(color="black", size=4, fontface="bold", box.padding=0.6,
                  point.padding=0.5, min.segment.length=0, max.overlaps=Inf) +
  scale_color_manual(values=cols_bio) + theme_classic() +
  ggtitle("P2. Colorectal (OSTA) Visium HD\nPseudo-bulk PCA (After)") +
  xlab(paste0("PC1 (", pb_res_after$pc1_var, "%)")) +
  ylab(paste0("PC2 (", pb_res_after$pc2_var, "%)")) +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5)) +
  scale_x_continuous(expand=expansion(mult=0.2)) +
  scale_y_continuous(expand=expansion(mult=0.2))
 
supp_plot <- (p_s1 + p_s2) /
             (p_u1 + p_u2) /
             (p_p1 + p_p2) +
             plot_layout(heights=c(1, 0.8, 0.8))
 
ggsave("Supp_Figure4_VisiumHD.png", supp_plot, width=12, height=15, dpi=300, bg="white")
 