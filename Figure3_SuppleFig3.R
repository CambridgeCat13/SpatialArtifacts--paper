library(SpatialExperiment)
library(SpatialArtifacts)
library(spatialLIBD)
library(scater)
library(scran)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(ggrepel)


sample_path <- "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/DLPFC_Br8325_ant_manual_alignment_all/outs/"
sample_id_name <- "DLPFC_Br8325"

spe_dlpfc <- read10xVisium(samples = sample_path, sample_id = sample_id_name, type = "HDF5", data = "filtered")
spe_dlpfc <- addPerCellQC(spe_dlpfc)
colData(spe_dlpfc)$sum_umi <- colSums(counts(spe_dlpfc))

spe_detected <- detectEdgeArtifacts(
  spe_dlpfc,
  platform = "visium",
  qc_metric = "sum",
  samples = "sample_id",
  batch_var = "sample_id",
  mad_threshold = 3,
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

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")
spe_ref <- spe[, spe$sample_id == "Br8325_ant"]
anno_vec <- colData(spe_ref)$bayesSpace_harmony_7
names(anno_vec) <- colnames(spe_ref)
colData(spe_classified)$cluster <- factor(
  anno_vec[match(colnames(spe_classified), names(anno_vec))]
)
spe_classified$cluster <- as.character(spe_classified$cluster)
spe_classified$cluster[is.na(spe_classified$cluster)] <- "Unannotated"
spe_classified$cluster <- factor(spe_classified$cluster)

plot_data_all <- cbind(
  as.data.frame(colData(spe_classified)),
  as.data.frame(spatialCoords(spe_classified)),
  as.data.frame(reducedDim(spe_classified, "UMAP"))
)
colnames(plot_data_all)[(ncol(plot_data_all)-1):ncol(plot_data_all)] <- c("UMAP1", "UMAP2")

plot_data_all$cluster <- as.character(plot_data_all$cluster)
plot_data_all$cluster[is.na(plot_data_all$cluster) | plot_data_all$cluster == "NA" | plot_data_all$cluster == "Unannotated"] <- "Unannotated"

all_levels <- sort(unique(plot_data_all$cluster))
cluster_levels <- c(setdiff(all_levels, "Unannotated"), "Unannotated")
plot_data_all$cluster <- factor(plot_data_all$cluster, levels = cluster_levels)

cls_col <- "edge_artifact_classification"
plot_data_all$class_factor <- factor(plot_data_all[[cls_col]], levels = c(
  "large_edge_artifact", "small_edge_artifact",
  "large_interior_artifact", "small_interior_artifact",
  "not_artifact"
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

n_total <- length(cluster_levels)
colors_colorful <- setNames(scales::hue_pal()(n_total), cluster_levels)
colors_colorful["Unannotated"] <- "#000000"

colors_highlight <- setNames(rep("#E0E0E0", n_total), cluster_levels)
colors_highlight["Unannotated"] <- "#D32F2F"

p_a <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=log10(sum_umi))) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_viridis(option="magma", name="log10(UMI)") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("A. Total UMI Distribution")

p_b <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=sum_umi < 1000)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("FALSE"="#F0F0F0", "TRUE"="black")) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("B. Standard QC (<1000)")

plot_data_all$draw_order <- ifelse(plot_data_all$class_factor == "Normal Spot", 1, 2)
p_c <- ggplot(plot_data_all[order(plot_data_all$draw_order),], aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=class_factor)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=artifact_colors) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("C. SpatialArtifacts Classification")

p_d <- ggplot(plot_data_all, aes(pxl_col_in_fullres, pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colors_colorful) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("D. Spatial Clusters (w/ Unannotated)")

p_e <- ggplot(plot_data_all[plot_data_all$class_factor == "Normal Spot",], aes(pxl_col_in_fullres, pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colors_colorful) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("E. Spatial Clusters (Filtered)")

p_f <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors, name="Artifact Type") +
  scale_y_continuous(labels=percent) +
  theme_classic() +
  ggtitle("F. Artifacts Explain Unannotated Spots") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom")

plot_data_all$hl_order <- ifelse(plot_data_all$cluster == "Unannotated", 2, 1)
p_g <- ggplot(plot_data_all[order(plot_data_all$hl_order),], aes(pxl_col_in_fullres, pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colors_highlight) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("G. Anatomy of Unannotated Spots")

p_h <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", color="black", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("H. UMI Distribution") + ylab("log10(UMI)") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_i <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", color="black", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("I. Gene Distribution") + ylab("log10(Genes)") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_j <- ggplot(plot_data_all[order(plot_data_all$hl_order),], aes(UMAP1, UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("J. UMAP (All Spots)")

p_k <- ggplot(plot_data_all[plot_data_all$class_factor == "Normal Spot",], aes(UMAP1, UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("K. UMAP (Filtered)")

p_l <- ggplot(plot_data_all, aes(UMAP1, UMAP2, color=sum_umi < 1000)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=c("FALSE"="#E0E0E0", "TRUE"="#D32F2F")) +
  theme_classic() +
  ggtitle("L. Global Threshold (<1000)")

final_plot_fig3 <- (p_a + p_b + p_c) / (p_d + p_e + p_f) / (p_g + p_h + p_i) / (p_j + p_k + p_l) +
  plot_layout(heights = c(1, 1, 0.8, 1))
ggsave("Figure3_DLPFC.png", plot=final_plot_fig3, width=15, height=20, dpi=300, bg="white")

# Supplementary Figure 3: Pseudo-bulk PCA

run_safe_pca <- function(spe_obj, title_text, color_scheme) {
  pb <- aggregateAcrossCells(spe_obj, ids = spe_obj$cluster)
  keep <- colSums(counts(pb)) > 100
  if (sum(keep) < 2) return(NULL)
  pb <- pb[, keep]

  sf <- librarySizeFactors(pb)
  sf[sf <= 0] <- min(sf[sf > 0])

  pb <- logNormCounts(pb, size.factors = sf)
  pb <- runPCA(pb, ncomponents = min(5, ncol(pb)-1))

  df <- as.data.frame(reducedDim(pb, "PCA"))
  df$cluster <- colnames(pb)
  var_exp <- attr(reducedDim(pb, "PCA"), "percentVar")
  pc1_v <- round(var_exp[1], 1)
  pc2_v <- round(var_exp[2], 1)

  ggplot(df, aes(PC1, PC2, color=cluster, label=cluster)) +
    geom_point(size=5, alpha=0.9) +
    geom_text_repel(fontface="bold", color="black", max.overlaps=Inf) +
    scale_color_manual(values=color_scheme) +
    theme_classic() +
    ggtitle(title_text) +
    xlab(paste0("PC1 (", pc1_v, "%)")) +
    ylab(paste0("PC2 (", pc2_v, "%)")) +
    theme(legend.position="none", plot.title=element_text(hjust=0.5, face="bold"))
}

p_supp_a <- run_safe_pca(
  spe_classified,
  "A. Before Artifact Removal (Including Unannotated)",
  colors_colorful
)

p_supp_b <- run_safe_pca(
  spe_classified[, spe_classified$edge_artifact_classification == "not_artifact" &
                   spe_classified$cluster != "Unannotated"],
  "B. After Artifact Removal",
  colors_colorful
)

supp_fig <- p_supp_a + p_supp_b +
  plot_annotation(
    title = "Supplementary Figure 3: Validation via Pseudo-bulk PCA",
    theme = theme(plot.title=element_text(hjust=0.5, face="bold", size=16))
  )

ggsave("Supp_Figure3_DLPFC.png", supp_fig, width=12, height=6, dpi=300, bg="white")