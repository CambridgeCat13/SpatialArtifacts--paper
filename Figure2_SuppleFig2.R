library(SpatialExperiment)
library(spatialLIBD)
library(scater)
library(scran)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(ggrepel)
library(SpatialArtifacts)

sample_path <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger-all/V11L05-335_C1/outs/"
sample_id_name <- "V11L05-335_C1"

spe <- read10xVisium(
  samples = sample_path,
  sample_id = sample_id_name,
  type = "HDF5",
  data = "filtered"
)
spe <- addPerCellQC(spe)
colData(spe)$sum <- colSums(counts(spe))

spe_detected <- detectEdgeArtifacts(
  spe,
  platform = "visium",
  qc_metric = "sum",
  samples = "sample_id",
  batch_var = "sample_id",
  mad_threshold = 2,
  edge_threshold = 0.5,
  name = "edge_artifact"
)

if (is.null(colData(spe_detected)$sum_umi)) {
  colData(spe_detected)$sum_umi <- colData(spe_detected)$sum
}

spe_classified <- classifyEdgeArtifacts(
  spe_detected,
  qc_metric = "sum_umi",
  min_spots = 20,
  name = "edge_artifact"
)

set.seed(123)
spe_classified <- logNormCounts(spe_classified)
spe_classified <- runPCA(spe_classified, ncomponents = 30)
spe_classified <- runUMAP(spe_classified)

g <- buildSNNGraph(spe_classified, k = 20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
spe_classified$cluster <- factor(clust)

plot_data_all <- cbind(
  as.data.frame(colData(spe_classified)),
  as.data.frame(spatialCoords(spe_classified)),
  as.data.frame(reducedDim(spe_classified, "UMAP"))
)
colnames(plot_data_all)[(ncol(plot_data_all)-1):ncol(plot_data_all)] <- c("UMAP1", "UMAP2")

cls_col <- "edge_artifact_classification"
plot_data_all$class_factor <- factor(plot_data_all[[cls_col]], levels = c(
  "large_edge_artifact", "small_edge_artifact",
  "large_interior_artifact", "small_interior_artifact",
  "not_artifact"
))
levels(plot_data_all$class_factor) <- c(
  "Large Edge", "Small Edge",
  "Large Interior", "Small Interior",
  "Normal Spot"
)

artifact_colors <- c(
  "Normal Spot"     = "lightgray",
  "Small Edge"      = "#FFB347",
  "Large Edge"      = "#FF4500",
  "Small Interior"  = "#00CED1",
  "Large Interior"  = "#0047AB"
)

real_clusters <- sort(unique(plot_data_all$cluster))
n_clusters <- length(real_clusters)

if (exists(".get_palette", where = asNamespace("scater"), mode = "function")) {
  colors_colorful <- scater:::.get_palette("tableau20")[1:n_clusters]
} else {
  colors_colorful <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                       "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
                       "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                       "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")[1:n_clusters]
}
names(colors_colorful) <- real_clusters

colors_highlight <- rep("#E0E0E0", n_clusters)
names(colors_highlight) <- real_clusters
if ("8" %in% names(colors_highlight)) colors_highlight["8"] <- "#FFB347"
if ("9" %in% names(colors_highlight)) colors_highlight["9"] <- "#FF4500"

p_a <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=log10(sum_umi))) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_viridis(option="magma", name="log10(UMI)") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("A. Total UMI Distribution") +
  theme(legend.position="right", plot.title=element_text(face="bold"))

std_threshold <- 1000
plot_data_all$std_qc_status <- ifelse(plot_data_all$sum_umi < std_threshold, "Fail", "Pass")
p_b <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=std_qc_status)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("Fail"="black", "Pass"="#F0F0F0"), name="QC") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("B. Standard QC (<1000)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

plot_data_all$draw_order <- ifelse(plot_data_all$class_factor == "Normal Spot", 1, 2)
plot_data_ordered <- plot_data_all[order(plot_data_all$draw_order), ]

p_c <- ggplot(plot_data_ordered, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=class_factor)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=artifact_colors) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("C. SpatialArtifacts Classification") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_d <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("D. Spatial Clusters (All)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

plot_data_filt <- plot_data_all[plot_data_all$class_factor == "Normal Spot", ]
p_e <- ggplot(plot_data_filt, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("E. Spatial Clusters (Filtered)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_f <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors, name="Artifact Type") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  ggtitle("F. Artifact Composition") +
  ylab("%") + xlab("Cluster ID") +
  theme(legend.position="bottom", legend.key.size=unit(0.4,"cm"), plot.title=element_text(face="bold"))

plot_data_all$highlight_order <- ifelse(plot_data_all$cluster %in% c("8", "9"), 2, 1)
plot_data_hl <- plot_data_all[order(plot_data_all$highlight_order), ]

p_g <- ggplot(plot_data_hl, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_highlight) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("G. Anatomy: Clusters 8 & 9") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_h <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("H. UMI Distribution") + ylab("log10(UMI)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_i <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("I. Gene Distribution") + ylab("log10(Genes)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_j <- ggplot(plot_data_hl, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("J. UMAP (All Spots)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_k <- ggplot(plot_data_filt, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("K. UMAP (Filtered)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

global_cutoff <- 1000
plot_data_all$is_global_low <- plot_data_all$sum_umi < global_cutoff
plot_data_global <- plot_data_all[order(plot_data_all$is_global_low), ]

p_l <- ggplot(plot_data_global, aes(x=UMAP1, y=UMAP2, color=is_global_low)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=c("FALSE"="#E0E0E0", "TRUE"="#D32F2F"),
                     labels=c("Pass", paste0("Fail (<", global_cutoff, ")"))) +
  theme_classic() +
  ggtitle(paste0("L. Global Threshold (UMI < ", global_cutoff, ")")) +
  theme(legend.position="bottom", legend.title=element_blank(), plot.title=element_text(face="bold"))

final_plot_fig2 <- (p_a + p_b + p_c) /
                   (p_d + p_e + p_f) /
                   (p_g + p_h + p_i) /
                   (p_j + p_k + p_l) +
                   plot_layout(heights = c(1, 1, 0.8, 1))

ggsave("Figure2_Hippo.png", plot=final_plot_fig2, width=15, height=20, dpi=300, bg="white")

# Supplementary Figure 2: Pseudo-bulk PCA

pb_all <- aggregateAcrossCells(spe_classified, ids = spe_classified$cluster)
keep_clusters <- colSums(counts(pb_all)) > 0
pb_all <- pb_all[, keep_clusters]
pb_all <- logNormCounts(pb_all, size.factors = librarySizeFactors(pb_all))
pb_all <- runPCA(pb_all, ncomponents = min(5, ncol(pb_all)-1))

pb_df_all <- as.data.frame(reducedDim(pb_all, "PCA"))
pb_df_all$cluster <- colnames(pb_all)

var_explained_all <- attr(reducedDim(pb_all, "PCA"), "percentVar")
pc1_var_all <- round(var_explained_all[1], 1)
pc2_var_all <- round(var_explained_all[2], 1)

p_supp_a <- ggplot(pb_df_all, aes(x=PC1, y=PC2, color=cluster, label=cluster)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(color="black", size=4, fontface="bold",
                  box.padding=0.6, point.padding=0.5,
                  min.segment.length=0, max.overlaps=Inf) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("A. Hippocampus (V11L05-335_C1)\nVisium - Before Artifact Removal") +
  xlab(paste0("PC1 (", pc1_var_all, "%)")) +
  ylab(paste0("PC2 (", pc2_var_all, "%)")) +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5)) +
  scale_x_continuous(expand=expansion(mult=0.2)) +
  scale_y_continuous(expand=expansion(mult=0.2))

spe_clean <- spe_classified[, spe_classified$edge_artifact_classification == "not_artifact"]

pb_clean <- aggregateAcrossCells(spe_clean, ids = spe_clean$cluster)
keep_clean <- colSums(counts(pb_clean)) > 0
pb_clean <- pb_clean[, keep_clean]

pb_clean <- logNormCounts(pb_clean, size.factors = librarySizeFactors(pb_clean))
pb_clean <- runPCA(pb_clean, ncomponents = min(5, ncol(pb_clean)-1))

pb_df_clean <- as.data.frame(reducedDim(pb_clean, "PCA"))
pb_df_clean$cluster <- colnames(pb_clean)

var_explained_clean <- attr(reducedDim(pb_clean, "PCA"), "percentVar")
pc1_var_clean <- round(var_explained_clean[1], 1)
pc2_var_clean <- round(var_explained_clean[2], 1)

p_supp_b <- ggplot(pb_df_clean, aes(x=PC1, y=PC2, color=cluster, label=cluster)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(color="black", size=4, fontface="bold",
                  box.padding=0.6, point.padding=0.5,
                  min.segment.length=0, max.overlaps=Inf) +
  scale_color_manual(values=colors_colorful) +
  theme_classic() +
  ggtitle("B. Hippocampus (V11L05-335_C1)\nVisium - After Artifact Removal") +
  xlab(paste0("PC1 (", pc1_var_clean, "%)")) +
  ylab(paste0("PC2 (", pc2_var_clean, "%)")) +
  theme(legend.position="none", plot.title=element_text(face="bold", size=12, hjust=0.5)) +
  scale_x_continuous(expand=expansion(mult=0.2)) +
  scale_y_continuous(expand=expansion(mult=0.2))

supp_fig <- p_supp_a + p_supp_b +
  plot_layout(ncol=2) +
  plot_annotation(
    title = "Supplementary Figure 2: Validation via Pseudo-bulk PCA",
    theme = theme(plot.title=element_text(size=18, face="bold", hjust=0.5))
  )

ggsave("Supp_Figure2_PseudoBulk.png", supp_fig, width=12, height=6, dpi=300, bg="white")