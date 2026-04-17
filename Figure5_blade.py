import os
import pandas as pd
import scanpy as sc
import time
from artifactsRemoval import Artifact_remove_HD


def prepare_standard_visium_as_hd(outs_dir):
    spatial_dir = os.path.join(outs_dir, "spatial")
    h5_path = os.path.join(outs_dir, "filtered_feature_bc_matrix.h5")
    if not os.path.exists(h5_path):
        h5_path = os.path.join(outs_dir, "filtered_feature_bc_matrix", "matrix.h5")

    csv_path = os.path.join(spatial_dir, "tissue_positions_list.csv")
    if not os.path.exists(csv_path):
        csv_path = os.path.join(spatial_dir, "tissue_positions.csv")

    target_parquet = os.path.join(spatial_dir, "tissue_positions_fixed.parquet")

    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path, header=None)
            first_val = str(df.iloc[0, 0])
            if "barcode" not in first_val.lower():
                df.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
            else:
                df = pd.read_csv(csv_path, header=0)
            df = df[df['in_tissue'] == 1]
            df.to_parquet(target_parquet, index=False)
            return h5_path, target_parquet
        except Exception as e:
            print(f"Error preparing Parquet: {e}")
            return None, None
    return None, None


hd_root = "/users/hjialihe/R_server/VisiumHD_Temp_Unzip/binned_outputs"
samples = [
    {
        "id": "1. Hippo (Standard)",
        "type": "Standard",
        "path": "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger-all/V11L05-335_C1/outs"
    },
    {
        "id": "2. DLPFC (Standard)",
        "type": "Standard",
        "path": "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger/DLPFC_Br8325_ant_manual_alignment_all/outs"
    },
    {
        "id": "3. OSTD Colorectal 16um (Visium HD)",
        "type": "HD",
        "h5_path": os.path.join(hd_root, "square_016um", "filtered_feature_bc_matrix.h5"),
        "pos_path": os.path.join(hd_root, "square_016um", "spatial/tissue_positions.parquet")
    }
]

output_dir = "/users/hjialihe/R_server/BLADE_Results_Dark"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for idx, sample in enumerate(samples):
    sample_id = sample["id"]

    if sample["type"] == "Standard":
        h5_path, pos_path = prepare_standard_visium_as_hd(sample["path"])
    else:
        h5_path, pos_path = sample["h5_path"], sample["pos_path"]

    k_neigh = 8 if sample["type"] == "Standard" else 6

    blade_obj = Artifact_remove_HD(
        data_dir=os.path.dirname(h5_path),
        positions_path=pos_path,
        h5_path=h5_path,
        neighbor_mode="k",
        k_neighbors=k_neigh
    )

    df_all = pd.read_parquet(pos_path)
    if 'barcode' in df_all.columns:
        df_all.set_index('barcode', inplace=True)
    if 'in_tissue' in df_all.columns:
        df_plot = df_all[df_all['in_tissue'] == 1].copy()
    else:
        df_plot = df_all.copy()

    radius = float(blade_obj._grid_radius) if hasattr(blade_obj, '_grid_radius') else 1.0
    blade_obj.remove_border(distance=radius)
    blade_obj.remove_edge(distance=2)
    try:
        blade_obj.remove_malfunction()
    except Exception as e_mal:
        print(f"WARNING: remove_malfunction failed: {e_mal}. Continuing without it.")

    temp_h5 = os.path.join(output_dir, f"temp_{idx}.h5")
    blade_obj.save_hdf5(temp_h5)

    adata_temp = sc.read_10x_h5(temp_h5)
    survivor_barcodes = set(adata_temp.obs_names)
    if os.path.exists(temp_h5):
        os.remove(temp_h5)

    df_plot['BLADE_classification'] = df_plot.index.map(
        lambda x: 'Normal' if x in survivor_barcodes else 'Artifact'
    )
    csv_path = os.path.join(output_dir, f"blade_labels_{idx}.csv")
    df_plot[['BLADE_classification']].to_csv(csv_path)