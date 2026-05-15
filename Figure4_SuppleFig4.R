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

p_g <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  scale_fill_manual(values=artifact_colors) +
  ggtitle("D. UMI Distribution") + ylab("log10(UMI)") + theme_vio

p_h <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  scale_fill_manual(values=artifact_colors) +
  ggtitle("E. Gene Distribution") + ylab("log10(Genes)") + theme_vio

target_cluster_id <- "1"
cluster1_data <- subset(plot_data_all, as.character(plot_data_all$cluster) == target_cluster_id)
other_data <- plot_data_all[as.character(plot_data_all$cluster) != target_cluster_id, ]

p_d <- ggplot() +
  geom_point(data=other_data, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres),
             color="grey95", size=hd_point_size, stroke=0) +
  geom_point(data=cluster1_data[order(ifelse(cluster1_data$class_factor == "Normal Spot", 1, 2)),],
             aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=class_factor),
             size=hd_point_size, stroke=0, alpha=1) +
  scale_color_manual(values=artifact_colors, name="Type") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle("F. Anatomy of Cluster 1") +
  my_theme + theme(legend.position="right")

p_e <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors) +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  ggtitle("G. Artifact Composition") + ylab("Percentage") + xlab("Cluster ID") +
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
  ggtitle(paste0("H. Global Threshold (UMI < ", global_cutoff, ")")) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        plot.title=element_text(size=16, face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        plot.margin=margin(2, 2, 2, 2))

final_plot_complete <- (p_a | p_b | p_c) /
                       (p_g | p_h) /
                       (p_d | p_e | p_f) +
                       plot_layout(heights=c(1.2, 0.9, 1.0)) &
                       theme(plot.margin=margin(1, 1, 1, 1))

ggsave("Figure4_VisiumHD_Colon.png", plot=final_plot_complete, 
       width=15, height=16, dpi=300, bg="white")

# Supplementary Figure 4
library(scater)
library(cowplot)
library(dplyr)

n_clust <- length(unique(plot_data_all$cluster))
cols_bio <- scales::hue_pal()(n_clust)
names(cols_bio) <- levels(plot_data_all$cluster)

df_clean <- plot_data_all[plot_data_all$class_factor == "Normal Spot", ]

p_s1 <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, 
                                   color=cluster)) +
  geom_point(size=hd_point_size, stroke=0, alpha=hd_alpha) +
  scale_color_manual(values=cols_bio, name="Cluster") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle("A. Spatial Clusters (All)") +
  theme_void() +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5),
        legend.position="bottom",
        legend.key.size=unit(0.3,"cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=8)) +
  guides(color=guide_legend(nrow=3, override.aes=list(size=3)))

p_s2 <- ggplot(df_clean, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, 
                               color=cluster)) +
  geom_point(size=hd_point_size, stroke=0, alpha=hd_alpha) +
  scale_color_manual(values=cols_bio, name="Cluster") +
  coord_fixed(xlim=x_limits, ylim=y_limits) + scale_y_reverse() +
  ggtitle("B. Spatial Clusters (Filtered)") +
  theme_void() +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5),
        legend.position="bottom",
        legend.key.size=unit(0.3,"cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=8)) +
  guides(color=guide_legend(nrow=3, override.aes=list(size=3)))

spe_classified$cluster <- factor(spe_classified$cluster)
spe_classified$edge_dryspot_hd_classification <- factor(
  spe_classified$edge_dryspot_hd_classification)

p_expl_before <- plotExplanatoryPCs(
  spe_classified,
  variables = c("sum_umi", "detected", "edge_dryspot_hd_classification"),
  npcs = 10
) + ggtitle("C. Explanatory PCs (Before Filtering)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

spe_clean_obj <- spe_classified[, spe_classified$edge_dryspot_hd_classification == "not_artifact"]
spe_clean_obj <- runPCA(spe_clean_obj, ncomponents = 30)

p_expl_after <- plotExplanatoryPCs(
  spe_clean_obj,
  variables = c("sum_umi", "detected"),
  npcs = 10
) + ggtitle("D. Explanatory PCs (After Filtering)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

plot_data_all$comparison_group <- case_when(
  plot_data_all$class_factor == "Large Edge"     ~ "Large Edge\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Small Edge"     ~ "Small Edge\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Large Interior" ~ "Large Interior\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Small Interior" ~ "Small Interior\n(SpatialArtifacts)",
  plot_data_all$sum_umi < 500 & plot_data_all$class_factor == "Normal Spot" ~ "False Positive\n(Standard QC only)",
  TRUE ~ "Normal Spot"
)

plot_data_all$comparison_group <- factor(plot_data_all$comparison_group, levels = c(
  "Large Edge\n(SpatialArtifacts)",
  "Small Edge\n(SpatialArtifacts)",
  "Large Interior\n(SpatialArtifacts)",
  "Small Interior\n(SpatialArtifacts)",
  "False Positive\n(Standard QC only)",
  "Normal Spot"
))

comparison_colors <- c(
  "Large Edge\n(SpatialArtifacts)"      = "#FF4500",
  "Small Edge\n(SpatialArtifacts)"      = "#FFB347",
  "Large Interior\n(SpatialArtifacts)"  = "#0047AB",
  "Small Interior\n(SpatialArtifacts)"  = "#00CED1",
  "False Positive\n(Standard QC only)"  = "#9370DB",
  "Normal Spot"                          = "lightgray"
)

p_comparison <- ggplot(plot_data_all,
                       aes(x = comparison_group, y = log10(sum_umi + 1),
                           fill = comparison_group)) +
  geom_violin(scale = "width", color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = comparison_colors) +
  theme_classic() +
  ggtitle("E. SpatialArtifacts vs Standard QC: UMI Distribution") +
  ylab("log10(UMI)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

epcam_id <- "ENSG00000119888"
acta2_id <- "ENSG00000107796"
cdh1_id  <- "ENSG00000039068"

genes_of_interest <- c("EPCAM", "ACTA2", "CDH1")
ensembl_ids <- c(epcam_id, acta2_id, cdh1_id)

for (i in seq_along(genes_of_interest)) {
  colData(spe_classified)[[genes_of_interest[i]]] <-
    as.numeric(logcounts(spe_classified)[ensembl_ids[i], ])
}

plot_data_genes <- as.data.frame(colData(spe_classified))

plot_list_hd <- lapply(genes_of_interest, function(gene) {
  ggplot(plot_data_genes,
         aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
             color = .data[[gene]])) +
    geom_point(size = hd_point_size, stroke = 0, alpha = hd_alpha) +
    scale_color_viridis_c(option = "magma", name = "log\nnorm") +
    coord_fixed(xlim = x_limits, ylim = y_limits) + scale_y_reverse() +
    theme_void() +
    theme(plot.title = element_text(face = "bold.italic", hjust = 0.5),
          legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
})
names(plot_list_hd) <- genes_of_interest

p_epcam <- plot_list_hd[["EPCAM"]] + ggtitle("F. EPCAM")
p_acta2 <- plot_list_hd[["ACTA2"]] + ggtitle("G. ACTA2")
p_cdh1  <- plot_list_hd[["CDH1"]]  + ggtitle("H. CDH1")        

row1 <- plot_grid(p_s1, p_s2, ncol=2)
row2 <- plot_grid(p_expl_before, p_expl_after, ncol=2)
row3 <- plot_grid(p_comparison, p_epcam, p_acta2, p_cdh1, ncol=4)

supp_fig4 <- plot_grid(
  row1, row2, row3,
  ncol=1,
  rel_heights=c(1.2, 1, 1)
)

ggsave("Supplementary_Figure4_VisiumHD.png", supp_fig4,
       width=18, height=16, dpi=300, bg="white")