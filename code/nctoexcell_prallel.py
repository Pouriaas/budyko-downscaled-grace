# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 22:48:06 2025

@author: pouria
"""
import numpy as np
# np.float = float
# np.round_ = np.round
import geopandas as gpd
import os
import pandas as pd
import warnings; warnings.filterwarnings(action='ignore')
import xarray as xr 

import regionmask
import pyogrio
from joblib import Parallel, delayed

warnings.filterwarnings(action='ignore')

# Load shapefile and initialize region mask
countries = gpd.read_file(r"/home/pouria/git/budyko-downscaled-grace/data/shapefile/Karkheh_karun.shp")

countries = countries.to_crs(4326)
indexes = list(range(len(countries)))
names=[str(i) for i in countries['Bas_Code'].values]
countries_mask_poly = regionmask.Regions(
    name='Name',
    numbers=indexes,
    names=names,
    abbrevs=names,
    outlines=list(countries.geometry.values)
)


ds = xr.open_mfdataset(r"/home/pouria/git/budyko-downscaled-grace/data/parametersinnc/era5.nc", decode_times=False)

# Main processing loop with parallel execution
mask = countries_mask_poly.mask(ds)
lat = mask.coords['lat'].values
lon = mask.coords['lon'].values
# Directory for output files
output_dir = r"/home/pouria/git/budyko-downscaled-grace/data/parametersinnc"

# Helper function to process each dataset
def process_dataset(ds_path, mask, lats, lons, ID_COUNTRY, c, country_code, lat1, lon1):
    ds = xr.open_mfdataset(ds_path, decode_times=False)
    sel_mask = mask.where(mask == ID_COUNTRY).values
    sel_mask[sel_mask == ID_COUNTRY] = 0
    id_lon = lons[np.where(~np.all(sel_mask, axis=0))]
    id_lat = lats[np.where(~np.all(sel_mask, axis=1))]
    
    if len(id_lon) > 0 and len(id_lat) > 0:
        out_sel1 = ds.sel(
            lat=slice(id_lat[0], id_lat[-1]),
            lon=slice(id_lon[0], id_lon[-1])
        ).compute().where(mask == ID_COUNTRY)
    else:
        out_sel1 = ds.sel(
            lat=lat1, lon=lon1, method="nearest"
        ).compute()
    
    var_names = list(out_sel1.variables.keys())
    for var in var_names:
        if var not in out_sel1.variables:
            print(f"Variable {var} not found in dataset {ds_path}. Skipping.")
            continue
        
        var_data = out_sel1[var]
        dims = var_data.dims  # Get dimensions of the variable
        
        if "lat" in dims and "lon" in dims:
            # For 3D variables with spatial dimensions
            mean_values = var_data.mean(dim=["lat", "lon"], skipna=True)
            c[var] = mean_values.values
        elif "lat" not in dims and "lon" not in dims:
            # For 1D variables or variables without spatial dimensions
            c[var] = var_data.values  # Use raw values
        else:
            print(f"Skipping variable {var}: Unsupported dimension structure.")
            continue

    return c
# Function to process a single country
def process_country(ID_COUNTRY, mask, lat, lon, countries, output_dir,ds):
    country_code = str(countries['Bas_Code'][ID_COUNTRY])
    print(f"Processing country: {country_code}")
    output_path = os.path.join(output_dir, f"{country_code}.xlsx")
    
    # Initialize dataframe
    centroid = countries.centroid.iloc[ID_COUNTRY]
    lon1, lat1 = centroid.x, centroid.y
    time = ds["time"].values
    # time_units = ds['time']
    # Extract units and origin from the dataset
    c = pd.DataFrame(time, columns=["time"])
    
    # Process each dataset
    datasets = [
        r"/home/pouria/git/budyko-downscaled-grace/data/parametersinnc/era5.nc",
    ]
    for dataset_path in datasets:
        c = process_dataset(dataset_path, mask, lat, lon, ID_COUNTRY, c, country_code, lat1, lon1)
    
    # Save to Excel
    c.to_excel(output_path, index=False)
    print(f"Finished processing country: {country_code}")



# Use Joblib to parallelize country processing
Parallel(n_jobs=11,verbose=5)(delayed(process_country)(ID_COUNTRY, mask, lat, lon, countries, output_dir,ds) for ID_COUNTRY in indexes)
