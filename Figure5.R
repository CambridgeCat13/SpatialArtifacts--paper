library(SpatialExperiment)
library(SpatialArtifacts)
library(SpotSweeper)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scuttle)
library(scater)
library(rhdf5)
library(arrow)

blade_dir <- "/users/hjialihe/R_server/BLADE_Results_Dark"

# VisiumHD(Colon)- MAD = 3

zip_file_path <- "/users/hjialihe/.cache/R/BiocFileCache/10049d62745472_VisiumHD_HumanColon_Oliveira.zip"
unzip_dir <- tempfile()
unzip(zip_file_path, exdir = unzip_dir)
h5_files <- list.files(unzip_dir, pattern = "filtered_feature_bc_matrix.h5", recursive = TRUE, full.names = TRUE)
h5_target <- h5_files[grep("square_016um", h5_files)]
if (length(h5_target) == 0) h5_target <- h5_files[1]

parquet_files <- list.files(unzip_dir, pattern = "tissue_positions.parquet", recursive = TRUE, full.names = TRUE)
parquet_target <- parquet_files[grep("square_016um", parquet_files)]
if (length(parquet_target) == 0) parquet_target <- parquet_files[1]
sce_temp <- read10xCounts(h5_target, col.names = TRUE)
xyz <- arrow::read_parquet(parquet_target)
match_idx <- match(colData(sce_temp)$Barcode, xyz$barcode)
xyz_sorted <- xyz[match_idx, ]

spe_hd <- SpatialExperiment(
  assays = list(counts = counts(sce_temp)),
  colData = colData(sce_temp),
  spatialCoords = as.matrix(xyz_sorted[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
)
colData(spe_hd)$array_row <- xyz_sorted$array_row
colData(spe_hd)$array_col <- xyz_sorted$array_col
colData(spe_hd)$in_tissue <- xyz_sorted$in_tissue
colData(spe_hd)$sample_id <- "VisiumHD_HumanColon_Oliveira"
spe_hd <- spe_hd[, which(colData(spe_hd)$in_tissue == 1)]
spe_hd <- addPerCellQC(spe_hd)
if ("sum" %in% colnames(colData(spe_hd))) colnames(colData(spe_hd))[colnames(colData(spe_hd)) == "sum"] <- "sum_umi"

spe_detected <- detectEdgeArtifacts(
  spe                  = spe_hd,
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
  keep_intermediate    = FALSE
)

spe_classified <- classifyEdgeArtifacts(
  spe       = spe_detected,
  qc_metric = "sum_umi",
  samples   = "sample_id",
  min_spots = 200,
  name      = "edge_dryspot_hd"
)
spe_hd_final <- spe_classified
spe_hd_final$res_ours <- ifelse(spe_hd_final$edge_dryspot_hd_classification == "not_artifact", "Normal", "Artifact")

#Hippocampus-MAD = 2

sample_path <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger-all/V11L05-335_C1/outs/"
sample_id_name <- "V11L05-335_C1"

spe_hippo <- read10xVisium(samples = sample_path, sample_id = sample_id_name, type = "HDF5", data = "filtered")
spe_hippo <- addPerCellQC(spe_hippo)
colData(spe_hippo)$sum_umi <- colSums(counts(spe_hippo))
spe_detected <- detectEdgeArtifacts(
  spe_hippo,
  platform      = "visium",
  qc_metric     = "sum",
  samples       = "sample_id",
  batch_var     = "sample_id",
  mad_threshold = 2,
  edge_threshold = 0.5,
  name          = "edge_artifact"
)
spe_classified <- classifyEdgeArtifacts(
  spe_detected,
  qc_metric = "sum_umi",
  min_spots = 20,
  name      = "edge_artifact"
)
spe_hippo_final <- spe_classified
spe_hippo_final$res_ours <- ifelse(spe_hippo_final$edge_artifact_classification == "not_artifact", "Normal", "Artifact")

#DLPFC-MAD=3
sample_path <- "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/DLPFC_Br8325_ant_manual_alignment_all/outs/"
sample_id_name <- "DLPFC_Br8325"
spe_dlpfc <- read10xVisium(samples = sample_path, sample_id = sample_id_name, type = "HDF5", data = "filtered")
spe_dlpfc <- addPerCellQC(spe_dlpfc)
colData(spe_dlpfc)$sum_umi <- colSums(counts(spe_dlpfc))

spe_detected <- detectEdgeArtifacts(
  spe_dlpfc,
  platform      = "visium",
  qc_metric     = "sum",
  samples       = "sample_id",
  batch_var     = "sample_id",
  mad_threshold = 3,
  edge_threshold = 0.5,
  name          = "edge_artifact"
)
spe_classified <- classifyEdgeArtifacts(
  spe_detected,
  qc_metric = "sum_umi",
  min_spots = 20,
  name      = "edge_artifact"
)
spe_dlpfc_final <- spe_classified
spe_dlpfc_final$res_ours <- ifelse(spe_dlpfc_final$edge_artifact_classification == "not_artifact", "Normal", "Artifact")


#SpotSweeper

clean_and_run_ss <- function(spe, k) {
  cols_to_remove <- c("sum_outliers", "local_outliers_sum", "z_score_sum")
  colData(spe) <- colData(spe)[, !(colnames(colData(spe)) %in% cols_to_remove)]
  if (is.null(colData(spe)$sum)) colData(spe)$sum <- colSums(counts(spe))
  spe <- localOutliers(spe, metric = "sum", direction = "lower", n_neighbors = k, log = TRUE, cutoff = 3)
  return(spe)
}
spe_hippo_final <- clean_and_run_ss(spe_hippo_final, k=36)
spe_dlpfc_final <- clean_and_run_ss(spe_dlpfc_final, k=36)
spe_hd_final    <- clean_and_run_ss(spe_hd_final, k=48)

#BLADE Results


load_blade_labels <- function(spe, csv_path) {
  blade_df <- read.csv(csv_path, row.names = 1)
  spe$blade_classification <- "Normal"
  common_barcodes <- intersect(colnames(spe), rownames(blade_df))
  artifact_barcodes <- rownames(blade_df)[blade_df$BLADE_classification == "Artifact"]
  spe$blade_classification[colnames(spe) %in% artifact_barcodes] <- "Artifact"
  return(spe)
}

spe_hippo_final <- load_blade_labels(spe_hippo_final, file.path(blade_dir, "blade_labels_0.csv"))
spe_dlpfc_final <- load_blade_labels(spe_dlpfc_final, file.path(blade_dir, "blade_labels_1.csv"))
spe_hd_final    <- load_blade_labels(spe_hd_final,    file.path(blade_dir, "blade_labels_2.csv"))

theme_dark_mode <- function() {
  theme_void() +
    theme(
      plot.background  = element_rect(fill="black", color=NA),
      panel.background = element_rect(fill="black", color=NA),
      text             = element_text(color="white", face="bold"),
      plot.title       = element_text(size=15, hjust=0.5, margin=margin(b=8)),
      legend.position  = "none"
    )
}

cols_neon <- c(
  "Artifact" = "#00FFFF",
  "Outlier"  = "#FF00FF",
  "BLADE"    = "#FFFF00",
  "Normal"   = "#A0A0A0"
)
plot_dark_panel <- function(spe, method_col, title, base_size=1.0) {
  df <- as.data.frame(colData(spe))
  coords <- as.data.frame(spatialCoords(spe))

  if ("pxl_col_in_fullres" %in% colnames(coords)) {
    df$x <- coords$pxl_col_in_fullres
    df$y <- coords$pxl_row_in_fullres
  } else {
    df$x <- coords$array_col
    df$y <- coords$array_row
  }

  if (method_col == "spatialdry") {
    df$Label <- ifelse(df$res_ours == "Artifact", "Artifact", "Normal")
  } else if (method_col == "spotsweeper") {
    df$Label <- ifelse(df$sum_outliers == TRUE, "Outlier", "Normal")
  } else if (method_col == "blade") {
    df$Label <- ifelse(df$blade_classification == "Artifact", "BLADE", "Normal")
  }

  df_bg <- df[df$Label == "Normal", ]
  df_fg <- df[df$Label != "Normal", ]

  ggplot() +
    geom_point(data=df_bg, aes(x=x, y=y),
               color=cols_neon["Normal"],
               size=base_size, stroke=0, alpha=0.7) +
    geom_point(data=df_fg, aes(x=x, y=y, color=Label),
               size=base_size * 2.5, stroke=0, alpha=1.0) +
    scale_color_manual(values=cols_neon) +
    scale_y_reverse() + coord_fixed() +
    ggtitle(title) +
    theme_dark_mode()
}

pt_size_std <- 1.5
pt_size_hd  <- 0.8
p1_sd <- plot_dark_panel(spe_hippo_final, "spatialdry",  "SpatialArtifacts\nHippocampus (V11L05-335_C1) Visium", pt_size_std)
p1_ss <- plot_dark_panel(spe_hippo_final, "spotsweeper", "SpotSweeper\nHippocampus (V11L05-335_C1) Visium", pt_size_std)
p1_bl <- plot_dark_panel(spe_hippo_final, "blade",       "BLADE\nHippocampus (V11L05-335_C1) Visium", pt_size_std)

p2_sd <- plot_dark_panel(spe_dlpfc_final, "spatialdry",  "\nDLPFC (Br8325_ant) Visium", pt_size_std)
p2_ss <- plot_dark_panel(spe_dlpfc_final, "spotsweeper", "\nDLPFC (Br8325_ant) Visium", pt_size_std)
p2_bl <- plot_dark_panel(spe_dlpfc_final, "blade",       "\nDLPFC (Br8325_ant) Visium", pt_size_std)

p3_sd <- plot_dark_panel(spe_hd_final, "spatialdry",  "\nColorectal (OSTA) Visium HD", pt_size_hd)
p3_ss <- plot_dark_panel(spe_hd_final, "spotsweeper", "\nColorectal (OSTA) Visium HD", pt_size_hd)
p3_bl <- plot_dark_panel(spe_hd_final, "blade",       "\nColorectal (OSTA) Visium HD", pt_size_hd)

final_fig5 <- (p1_sd | p1_ss | p1_bl) /
              (p2_sd | p2_ss | p2_bl) /
              (p3_sd | p3_ss | p3_bl) +
              plot_annotation(
                title = "Figure 5: Benchmarking Validation",
                theme = theme(
                  plot.background = element_rect(fill="black", color=NA),
                  plot.title = element_text(color="white", size=20, face="bold", hjust=0.5)
                )
              )

ggsave("Figure5_Benchmark.png", final_fig5, width=18, height=18, dpi=300, bg="black")