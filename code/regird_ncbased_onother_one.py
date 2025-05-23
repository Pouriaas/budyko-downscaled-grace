#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 11:16:30 2025

@author: pouria
"""

import xarray as xr
import xesmf as xe
import os

# Open the source and target datasets
target_ds = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/input_data/GRACE_monthly01.nc")
source_ds = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/input_data/era5_monthly.nc")

# Show variables and coordinate names to adjust if necessary
print("Source dataset variables:", source_ds.data_vars)
print("Target dataset variables:", target_ds.data_vars)
print("Source coords:", source_ds.coords)
print("Target coords:", target_ds.coords)

# Optional: Rename coordinates if needed (e.g., 'lon' instead of 'longitude')
# source_ds = source_ds.rename({'longitude': 'lon', 'latitude': 'lat'})
# target_ds = target_ds.rename({'longitude': 'lon', 'latitude': 'lat'})

# Build regridder with conservative method
regridder = xe.Regridder(source_ds, target_ds, method='conservative')

# Regrid a specific variable (e.g., 'precip')
regridded_var = regridder(source_ds)

# Save the regridded data
# out_put=os.makedirs (r"/home/pouria/git/budyko-downscaled-grace/data/out_put")
out_put=r"/home/pouria/git/budyko-downscaled-grace/data/out_put"
path=os.path.join(out_put,"regrided_era5_based_on_grace.nc")
regridded_var.to_netcdf(path)