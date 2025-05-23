# -*- coding: utf-8 -*-
"""
Crop ERA5 NetCDF data based on shapefile polygons and save each as a separate NetCDF file.
Created on Jan 5, 2025
@author: pouria
"""

import os
import geopandas as gpd
import xarray as xr
import rioxarray
from shapely.geometry import mapping
from joblib import Parallel, delayed
import warnings

warnings.filterwarnings("ignore")

# === Input paths ===
shapefile_path = "/home/pouria/git/budyko-downscaled-grace/data/shapefile/Karkheh_karun.shp"
netcdf_path = "/home/pouria/git/budyko-downscaled-grace/data/input_data/era5_monthly.nc"
output_dir = "/home/pouria/git/budyko-downscaled-grace/data/input_data"

os.makedirs(output_dir, exist_ok=True)

# === Load shapefile ===
regions = gpd.read_file(shapefile_path).to_crs(4326)  # Ensure CRS matches NetCDF
region_ids = regions["Bas_Code"].astype(str).tolist()

# === Load NetCDF with rioxarray ===
ds = xr.open_dataset(netcdf_path, decode_times=False)
ds = ds.rio.write_crs("EPSG:4326", inplace=True)  # Set CRS if not present
ds = ds.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")

# === Function to crop a single region ===
def crop_region(i):
    region = regions.iloc[i]
    geom = [mapping(region.geometry)]
    region_id = region_ids[i]
    
    print(f"Processing region: {region_id}")


    clipped = ds.rio.clip(geom, all_touched=True, drop=True, invert=False)
    output_path = os.path.join(output_dir, "era5_cut_by_basin.nc")
    clipped.to_netcdf(output_path)
    print(f"Saved: {output_path}")


# === Parallel crop processing ===
Parallel(n_jobs=8, verbose=5)(delayed(crop_region)(i) for i in range(len(regions)))