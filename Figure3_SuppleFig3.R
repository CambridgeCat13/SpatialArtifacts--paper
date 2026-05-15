library(SpatialExperiment)
library(SpatialArtifacts)
library(spatialLIBD)
library(scater)
library(scran)
library(ggplot2)
library(cowplot)
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
plot_data_all$cluster[is.na(plot_data_all$cluster) | plot_data_all$cluster == "NA"] <- "Unannotated"

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

p_a <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                                  color=log10(sum_umi))) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_viridis(option="magma", name="log10(UMI)") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("A. Total UMI Distribution") +
  theme(legend.position=c(0.15, 0.25),
        legend.background=element_rect(fill=alpha("white", 0.7)),
        legend.key.size=unit(0.3,"cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        plot.title=element_text(face="bold"))

p_b <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                                  color=sum_umi < 1000)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=c("FALSE"="#F0F0F0", "TRUE"="black")) +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("B. Standard QC (<1000)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

plot_data_all$draw_order <- ifelse(plot_data_all$class_factor == "Normal Spot", 1, 2)
plot_data_ordered <- plot_data_all[order(plot_data_all$draw_order), ]

p_c <- ggplot(plot_data_ordered, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                                      color=class_factor)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=artifact_colors, name="Classification") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("C. SpatialArtifacts Classification") +
  theme(legend.position="right", plot.title=element_text(face="bold")) +
  guides(color=guide_legend(override.aes=list(size=3)))

p_d <- ggplot(plot_data_all, aes(x=class_factor, y=log10(sum_umi+1),
                                  fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("D. UMI Distribution") + ylab("log10(UMI)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold"))

p_e <- ggplot(plot_data_all, aes(x=class_factor, y=log10(detected+1),
                                  fill=class_factor)) +
  geom_violin(scale="width", color=NA, alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
  scale_fill_manual(values=artifact_colors) +
  theme_classic() +
  ggtitle("E. Gene Distribution") + ylab("log10(Genes)") +
  theme(legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold"))

p_f <- ggplot(plot_data_all, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                                  color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful, name="Cluster") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("F. Spatial Clusters (All)") +
  theme(legend.position="bottom", plot.title=element_text(face="bold"),
        legend.key.size=unit(0.4,"cm")) +
  guides(color=guide_legend(nrow=2, override.aes=list(size=3)))

plot_data_filt <- plot_data_all[plot_data_all$class_factor == "Normal Spot", ]
p_g <- ggplot(plot_data_filt, aes(x=pxl_col_in_fullres, y=pxl_row_in_fullres,
                                   color=cluster)) +
  geom_point(size=0.8, alpha=0.9) +
  scale_color_manual(values=colors_colorful, name="Cluster") +
  theme_void() + coord_fixed() + scale_y_reverse() +
  ggtitle("G. Spatial Clusters (Filtered)") +
  theme(legend.position="bottom", plot.title=element_text(face="bold"),
        legend.key.size=unit(0.4,"cm")) +
  guides(color=guide_legend(nrow=2, override.aes=list(size=3)))

p_h <- ggplot(plot_data_all, aes(x=cluster, fill=class_factor)) +
  geom_bar(position="fill", width=0.8) +
  scale_fill_manual(values=artifact_colors, name="Artifact Type") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic() +
  ggtitle("H. Artifact Composition") +
  ylab("%") + xlab("Cluster") +
  theme(legend.position="bottom", legend.key.size=unit(0.4,"cm"),
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(face="bold"))

plot_data_all$hl_order <- ifelse(plot_data_all$cluster == "Unannotated", 2, 1)
plot_data_hl <- plot_data_all[order(plot_data_all$hl_order), ]

p_i <- ggplot(plot_data_hl, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("I. UMAP (All Spots)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

p_j <- ggplot(plot_data_filt, aes(x=UMAP1, y=UMAP2, color=cluster)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=colors_highlight) +
  theme_classic() +
  ggtitle("J. UMAP (Filtered)") +
  theme(legend.position="none", plot.title=element_text(face="bold"))

plot_data_all$is_global_low <- plot_data_all$sum_umi < 1000
plot_data_global <- plot_data_all[order(plot_data_all$is_global_low), ]

p_k <- ggplot(plot_data_global, aes(x=UMAP1, y=UMAP2, color=is_global_low)) +
  geom_point(size=0.5, alpha=0.6) +
  scale_color_manual(values=c("FALSE"="#E0E0E0", "TRUE"="#D32F2F"),
                     labels=c("Pass", "Fail (<1000)")) +
  theme_classic() +
  ggtitle("K. Global Threshold (UMI < 1000)") +
  theme(legend.position="bottom", legend.title=element_blank(),
        plot.title=element_text(face="bold"))

top_row  <- plot_grid(p_a, p_b, p_c, ncol=3, rel_widths=c(1,1,1))
mid_row1 <- plot_grid(p_d, p_e, ncol=2)
mid_row2 <- plot_grid(p_f, p_g, p_h, ncol=3, rel_widths=c(1.2,1.2,0.8))
bot_row  <- plot_grid(p_i, p_j, p_k, ncol=3)

final_plot_fig3 <- plot_grid(
  top_row, mid_row1, mid_row2, bot_row,
  ncol=1,
  rel_heights=c(1, 1, 1.2, 1)
)

ggsave("Figure3_DLPFC.png", plot=final_plot_fig3, width=15, height=24,
       dpi=300, bg="white")
# Supplementary Figure 3
library(scater)

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
  plot_data_all$class_factor == "Large Edge"     ~ "Large Edge\n(SpatialArtifacts)",
  plot_data_all$class_factor == "Small Edge"     ~ "Small Edge\n(SpatialArtifacts)",
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

snap25_id <- "ENSG00000132639"
gfap_id   <- "ENSG00000131095"
mobp_id   <- "ENSG00000168314"

genes_of_interest <- c("SNAP25", "GFAP", "MOBP")
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

top_row_supp <- plot_grid(p_expl_before, p_expl_after, ncol = 2)
bot_row_supp <- plot_grid(p_comparison, p_snap25, p_gfap, p_mobp, ncol = 4)

supp_fig3 <- plot_grid(
  top_row_supp,
  bot_row_supp,
  ncol = 1,
  rel_heights = c(1, 1)
)

ggsave("Supplementary_Figure3_DLPFC.png", supp_fig3,
       width = 18, height = 10, dpi = 300, bg = "white")