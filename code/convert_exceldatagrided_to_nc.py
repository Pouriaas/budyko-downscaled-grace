import pandas as pd
import numpy as np
import xarray as xr

# Step 1: Load Excel file
file_path = "/home/pouria/git/budyko-downscaled-grace/data/budko_parameters/0.1-nc4.xlsx"
df = pd.read_excel(file_path)



# Step 2: Extract coordinates and time
lons = df.iloc[:, 0].values
lats = df.iloc[:, 1].values

# Convert date columns to datetime
date_columns = df.columns[2:]
dates = pd.to_datetime(date_columns, format="%m/%d/%Y")

# Step 3: Unique sorted lat/lon values
unique_coords = df.iloc[:, :2].drop_duplicates()
lon_values = np.sort(unique_coords.iloc[:, 0].unique())
lat_values = np.sort(unique_coords.iloc[:, 1].unique())

# Step 4: Initialize data array (time x lat x lon)
data = np.full((len(dates), len(lat_values), len(lon_values)), np.nan)

# Step 5: Fill the data array
for _, row in df.iterrows():
    lon, lat = row.iloc[0], row.iloc[1]
    for i, date in enumerate(date_columns):
        val = row[date]
        if pd.notna(val):
            lat_idx = np.where(lat_values == lat)[0][0]
            lon_idx = np.where(lon_values == lon)[0][0]
            data[i, lat_idx, lon_idx] = val

# Step 6: Create xarray Dataset
ds = xr.Dataset(
    {"value": (["time", "lat", "lon"], data)},
    coords={"time": dates, "lat": lat_values, "lon": lon_values}
)

# Step 7: Save to NetCDF
output_path = "/home/pouria/git/budyko-downscaled-grace/data/input_data/GRACE_monthly01.nc"
ds.to_netcdf(output_path)