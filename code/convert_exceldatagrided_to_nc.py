import pandas as pd
import numpy as np
import xarray as xr

# Step 1: Load Excel file
df = pd.read_csv("/home/pouria/git/budyko-downscaled-grace/data/budko_parameters/GRACE_Era5/Era_land_Pre_aev_pot/param.csv")

df=df.drop(columns=df.columns[0])
# Step 2: Extract lon, lat, and years
lons = df.iloc[:, 0].values
lats = df.iloc[:, 1].values
years = df.columns[2:].astype(int)  # Assuming all remaining columns are years

# Step 3: Create mesh of lon-lat pairs
unique_coords = df.iloc[:, :2].drop_duplicates()
lon_values = np.sort(unique_coords.iloc[:, 0].unique())
lat_values = np.sort(unique_coords.iloc[:, 1].unique())

# Step 4: Initialize data array
data = np.full((len(years), len(lat_values), len(lon_values)), np.nan)

# Fill the data array
for i, year in enumerate(years):
    for _, row in df.iterrows():
        lon, lat = row.iloc[0], row.iloc[1]
        val = row[str(year)]
        lat_idx = np.where(lat_values == lat)[0][0]
        lon_idx = np.where(lon_values == lon)[0][0]
        data[i, lat_idx, lon_idx] = val

# Step 5: Create xarray Dataset
ds = xr.Dataset(
    {
        "value": (["time", "lat", "lon"], data)
    },
    coords={
        "time": years,
        "lat": lat_values,
        "lon": lon_values
    }
)

# Step 6: Save to NetCDF
ds.to_netcdf("/home/pouria/git/budyko-downscaled-grace/data/parametersinnc/era5.nc")
