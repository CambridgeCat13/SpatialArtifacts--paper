# SpatialArtifacts-paper

Code to reproduce figures and analyses for the SpatialArtifacts manuscript.

## Overview

This repository contains the R and Python scripts used to generate all computational figures in the manuscript. Figure 1 and Supplementary Figure 1 are schematic diagrams and are not included here.

## Requirements

### R packages
- SpatialArtifacts
- SpatialExperiment
- SpotSweeper
- scater, scran, scuttle
- ggplot2, patchwork, scales, viridis, ggrepel
- spatialLIBD
- arrow, rhdf5

### Python packages (for Figure 5 BLADE benchmarking only)
- scanpy
- pandas, numpy
- artifactsRemoval (BLADE)

Install BLADE following the instructions at: https://github.com/KummerfeldLab/BLADE

## Figure Scripts

| Script | Figure |
|--------|--------|
| `Figure2_SuppleFig2.R` | Figure 2 + Supplementary Figure 2 |
| `Figure3_SuppleFig3.R` | Figure 3 + Supplementary Figure 3 |
| `Figure4_SuppleFig4.R` | Figure 4 + Supplementary Figure 4 |
| `Figure5.R` | Figure 5 |
| `SuppleFig5_6.R` | Supplementary Figure 5 + Supplementary Figure 6 |
| `Figure5_blade.py` | BLADE results for Figure 5 |

## Running Order

Most figures can be run independently. The exception is Figure 5, which requires BLADE results to be generated first:

**Step 1:** Activate the BLADE conda environment and run the Python script:
```bash
conda activate your_blade_envirment
python Figure5_blade.py
```

**Step 2:** Run the R script:
```bash
Rscript Figure5.R
```

All other figures can be run directly with:
```bash
Rscript Figure2_SuppleFig2.R
Rscript Figure3_SuppleFig3.R
Rscript Figure4_SuppleFig4.R
Rscript SuppleFig5_6.R
```

## Data

Data used in this manuscript are publicly available. See the Data Availability section of the manuscript for details.