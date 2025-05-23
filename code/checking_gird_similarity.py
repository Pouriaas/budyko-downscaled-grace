#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 13:22:51 2025

@author: pouria
"""

import xarray as xr
import numpy as np

# Open both datasets
ds1 = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/out_put/regrided_era5_based_on_grace.nc")
ds2 = xr.open_dataset("/home/pouria/git/budyko-downscaled-grace/data/input_data/GRACE_monthly01.nc")

# Check coordinate names
# print("ds1 coords:", ds1.coords)
# print("ds2 coords:", ds2.coords)
# hkjnj

# Use your dataset's actual coordinate names
lat_equal = np.array_equal(ds1['lat'], ds2['lat'])
lon_equal = np.array_equal(ds1['lon'], ds2['lon'])

if lat_equal and lon_equal:
    print("✅ Datasets are on the same grid.")
else:
    print("❌ Datasets are NOT on the same grid.")