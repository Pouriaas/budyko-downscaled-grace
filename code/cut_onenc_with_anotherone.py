#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 18:17:50 2025

@author: pouria
"""

import xarray as xr

# Load GRACE file and extract mask (should be 2D: lat/lon)
grace_ds = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/input_data/GRACE_monthly01.nc")
grace_mask = ~grace_ds["value"].isel(time=0).isnull()  # 2D mask: True inside basin

# Load target dataset
target_ds = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/out_put/regrided_era5_based_on_grace.nc")

# Automatically apply the 2D mask across all time steps
clipped = target_ds.where(grace_mask)

# Save result
clipped.to_netcdf("/home/pouria/git/budyko-downscaled-grace/data/out_put/clipped_target.nc")