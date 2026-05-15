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
library(cowplot)

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

colors_colorful <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                     "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
                     "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                     "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")[1:n_clusters]
names(colors_colorful) <- real_clusters

# Cluster 8 and 9 get hot distinct colors
colors_colorful["8"] <- "#C2185B"   # hot pink/magenta
colors_colorful["9"] <- "#FF0000"   # deep amber

colors_highlight <- rep("#E0E0E0", n_clusters)
names(colors_highlight) <- real_clusters
colors_highlight["8"] <- "#C2185B"
colors_highlight["9"] <- "#FF0000"

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
  scale_color_manual(values=artifact_colors, name="Classification") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("C. SpatialArtifacts Classification") +
  theme(legend.position="right", plot.title=element_text(face="bold")) +
  guides(color=guide_legend(override.aes=list(size=3)))

p_h <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("D. UMI Distribution") + ylab("log10(UMI)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_i <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1), fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("E. Gene Distribution") + ylab("log10(Genes)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(face="bold"))

p_d <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful, name="Cluster") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("F. Spatial Clusters (All)") +
  theme(legend.position="bottom", plot.title=element_text(face="bold"),
        legend.key.size=unit(0.4,"cm")) +
  guides(color=guide_legend(nrow=2, override.aes=list(size=3)))

plot_data_filt <- plot_data_all[plot_data_all$class_factor == "Normal Spot", ]
p_e <- ggplot(plot_data_filt, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres, color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful, name="Cluster") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("G. Spatial Clusters (Filtered)") +
  theme(legend.position="bottom", plot.title=element_text(face="bold"),
        legend.key.size=unit(0.4,"cm")) +
  guides(color=guide_legend(nrow=2, override.aes=list(size=3)))

p_f <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors, name="Artifact Type") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  ggtitle("H. Artifact Composition") +
  ylab("%") + xlab("Cluster ID") +
  theme(legend.position="bottom", legend.key.size=unit(0.4,"cm"),
        plot.title=element_text(face="bold"))

plot_data_hl <- plot_data_all
plot_data_hl$highlight_order <- ifelse(plot_data_hl$cluster %in% c("8","9"), 2, 1)
plot_data_hl <- plot_data_hl[order(plot_data_hl$highlight_order), ]

p_j <- ggplot(plot_data_hl, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("I. UMAP (All Spots)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_k <- ggplot(plot_data_filt, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("J. UMAP (Filtered)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

global_cutoff <- 1000
plot_data_all$is_global_low <- plot_data_all$sum_umi < global_cutoff
plot_data_global <- plot_data_all[order(plot_data_all$is_global_low), ]

p_l <- ggplot(plot_data_global, aes(x=UMAP1, y=UMAP2, color=is_global_low)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=c("FALSE"="#E0E0E0", "TRUE"="#D32F2F"),
                     labels=c("Pass", paste0("Fail (<", global_cutoff, ")"))) +
  theme_classic() +
  ggtitle(paste0("K. Global Threshold (UMI < ", global_cutoff, ")")) +
  theme(legend.position="bottom", legend.title=element_blank(),
        plot.title=element_text(face="bold"))

top_row  <- plot_grid(p_a, p_b, p_c, ncol=3, rel_widths=c(1,1,1))
mid_row1 <- plot_grid(p_h, p_i, ncol=2, rel_widths=c(1,1))
mid_row2 <- plot_grid(p_d, p_e, p_f, ncol=3, rel_widths=c(1.2,1.2,0.8))
bot_row  <- plot_grid(p_j, p_k, p_l, ncol=3, rel_widths=c(1,1,1))

final_plot_fig2 <- plot_grid(
  top_row, mid_row1, mid_row2, bot_row,
  ncol=1,
  rel_heights=c(1, 1, 1.2, 1)
)

ggsave("Figure2_Hippo.png", plot=final_plot_fig2, width=15, height=24, dpi=300, bg="white")

#Supplementary Figure 2:
library(scater)
library(cowplot)
library(dplyr)

spe_classified$cluster <- factor(spe_classified$cluster)
spe_classified$edge_artifact_classification <- factor(spe_classified$edge_artifact_classification)

p_expl_before <- plotExplanatoryPCs(
  spe_classified,
  variables = c("sum_umi", "detected", "edge_artifact_classification"),
  npcs = 10
) + ggtitle("A. Explanatory PCs (Before Filtering)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

spe_clean <- spe_classified[, spe_classified$edge_artifact_classification == "not_artifact"]
spe_clean <- runPCA(spe_clean, ncomponents = 30)

p_expl_after <- plotExplanatoryPCs(
  spe_clean,
  variables = c("sum_umi", "detected"),
  npcs = 10
) + ggtitle("B. Explanatory PCs (After Filtering)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

plot_data_all$comparison_group <- case_when(
  plot_data_all$class_factor == "Large Edge" ~ "Large Edge\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Small Edge" ~ "Small Edge\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Small Interior" ~ "Small Interior\n(SpatialArtifacts)",
  plot_data_all$sum_umi < 1000 & plot_data_all$class_factor == "Normal Spot" ~ "False Positive\n(Standard QC only)",
  TRUE ~ "Normal Spot"
)

plot_data_all$comparison_group <- factor(plot_data_all$comparison_group, levels = c(
  "Large Edge\n(SpatialArtifacts)",
  "Small Edge\n(SpatialArtifacts)",
  "Small Interior\n(SpatialArtifacts)",
  "False Positive\n(Standard QC only)",
  "Normal Spot"
))

comparison_colors <- c(
  "Large Edge\n(SpatialArtifacts)"      = "#FF4500",
  "Small Edge\n(SpatialArtifacts)"      = "#FFB347",
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
  ggtitle("C. SpatialArtifacts vs Standard QC: UMI Distribution") +
  ylab("log10(UMI)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

genes_of_interest <- c("SNAP25", "GFAP", "MOBP")
snap25_id <- rownames(spe_classified)[rowData(spe_classified)$symbol == "SNAP25"]
gfap_id   <- rownames(spe_classified)[rowData(spe_classified)$symbol == "GFAP"]
mobp_id   <- rownames(spe_classified)[rowData(spe_classified)$symbol == "MOBP"]
ensembl_ids <- c(snap25_id, gfap_id, mobp_id)

for (i in seq_along(genes_of_interest)) {
  colData(spe_classified)[[genes_of_interest[i]]] <-
    as.numeric(logcounts(spe_classified)[ensembl_ids[i], ])
}

plot_data_genes <- cbind(
  as.data.frame(colData(spe_classified)),
  as.data.frame(spatialCoords(spe_classified))
)

plot_list <- lapply(genes_of_interest, function(gene) {
  ggplot(plot_data_genes,
         aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
             color = .data[[gene]])) +
    geom_point(size = 0.8, alpha = 0.9) +
    scale_color_viridis(option = "magma", name = "log\nnorm") +
    theme_void() + coord_fixed() + scale_y_reverse() +
    theme(plot.title = element_text(face = "bold.italic", hjust = 0.5),
          legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
})
names(plot_list) <- genes_of_interest


p_snap25 <- plot_list[["SNAP25"]] + ggtitle("D. SNAP25")
p_gfap   <- plot_list[["GFAP"]]   + ggtitle("E. GFAP")
p_mobp   <- plot_list[["MOBP"]]   + ggtitle("F. MOBP")

p_300 <- ggplot(plot_data_all,
                aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                    color=sum_umi < 300)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("FALSE"="#F0F0F0", "TRUE"="black")) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("G. Standard QC (UMI < 300)") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        legend.position="none")

p_500 <- ggplot(plot_data_all,
                aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                    color=sum_umi < 500)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("FALSE"="#F0F0F0", "TRUE"="black")) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("H. Standard QC (UMI < 500)") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        legend.position="none")

p_1000 <- ggplot(plot_data_all,
                 aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                     color=sum_umi < 1000)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("FALSE"="#F0F0F0", "TRUE"="black")) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("I. Standard QC (UMI < 1000)") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        legend.position="none")

threshold_row <- plot_grid(p_300, p_500, p_1000, ncol=3)

top_row_supp <- plot_grid(p_expl_before, p_expl_after, ncol=2)
bot_row_supp <- plot_grid(p_comparison, p_snap25, p_gfap, p_mobp, ncol=4)

supp_fig2 <- plot_grid(
  top_row_supp,
  bot_row_supp,
  threshold_row,
  ncol=1,
  rel_heights=c(1, 1, 0.8)
)

ggsave("Supp_Figure2.png", supp_fig2, width=18, height=14,
       dpi=300, bg="white")