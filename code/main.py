import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import spotpy
from spotpy.parameter import Uniform
from hydroeval import kge as hydro_kge

# Placeholder: implement your Budyko-based transfer function here
# This should replicate the logic of the R function Fu() from with_corroloatory.R
def Fu(alpha, v, c, d, e, f, NDVI, lon, lat, alt, curve, rain, pot):
    # TODO: translate the R logic into Python
    # Example return: simulated evapotranspiration
    raise NotImplementedError("Please implement the Budyko function logic")

# Load comparison dataframe exported from R (e.g., CSV)
data = pd.read_csv('data/comparison.xlsx')

# Optionally remove specified stations by code
omitted_codes = [
    "35-009","21-223","21-105","41-907","41-111","44-035","41-083",
    "21-261","16-003","12-027","21-131","21-273","34-011","44-105",
    "41-039","18-061","22-045","41-149","18-089","41-143","41-117",
    "21-275","33-015","41-165","15-015","16-051","16-021","12-043",
    "12-085","16-023","13-005","34-003","41-115","46-011","41-159",
    "47-029","31-013","41-099"
]
data = data[~data['Code'].isin(omitted_codes)].reset_index(drop=True)

# Define variables to normalize
features = ['NDVI', 'DEM', 'curve number', 'lon', 'lat']
scaler = MinMaxScaler()
data[features] = scaler.fit_transform(data[features])

# Split into train/test
train_frac = 0.7
train = data.sample(frac=train_frac, random_state=42)
test = data.drop(train.index).reset_index(drop=True)
train = train.reset_index(drop=True)

# Prepare observed series for train
rain_train = train['Rain_prod_vol'] - train['delta_storage']
pot_evap_train = train['Potential_Eva_prod_vol']
obs_evap_train = train['Eva_prod_vol']
obs_strflow_train = train['hydr_vol']

# SPOTPY setup
class BudykoSPOTPY:
    def __init__(self, data):
        self.data = data
        # Parameter bounds as in R: six decision variables between -20 and 20
        self.params = [Uniform('dec{}'.format(i+1), low=-20, high=20) for i in range(6)]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # Unpack decision variables
        dec = vector
        sim_evap = []
        sim_strflow = []
        # Loop through samples
        for i, row in self.data.iterrows():
            rain = row['Rain_prod_vol'] - row['delta_storage']
            pot = row['Potential_Eva_prod_vol']
            ndvi, lon, lat, alt, curve = (
                row['NDVI'], row['lon'], row['lat'], row['DEM'], row['curve number']
            )
            evap = Fu(dec[0], dec[1], dec[2], dec[3], dec[4], dec[5],
                      ndvi, lon, lat, alt, curve, rain, pot)
            sim_evap.append(evap)
            sim_strflow.append(rain - evap)
        return np.array(sim_evap), np.array(sim_strflow)

    def evaluation(self):
        ev_train = self.data['Eva_prod_vol']
        sf_train = self.data['hydr_vol']
        return ev_train.values, sf_train.values

    def objectivefunction(self, simulation, evaluation):
        sim_evap, sim_strflow = simulation
        obs_evap, obs_strflow = evaluation
        # Compute KGE for each
        kge_str = hydro_kge(sim_strflow, obs_strflow)
        kge_ev = hydro_kge(sim_evap, obs_evap)
        # Composite objective as in R: (KGE_str-1)*3 + (KGE_ev-1)
        return (kge_str - 1) * 3 + (kge_ev - 1)

# Initialize setup with training data
spot_setup = BudykoSPOTPY(train)

# Run SCE-UA
sampler = spotpy.algorithms.sceua(spot_setup,
                                 dbname='sceua_budyko',
                                 dbformat='csv',
                                 alt_objfun=True,
                                 sim_timeout=60)
results = sampler.sample(2000)  # adjust number of function evaluations

# Extract best parameters
best_idx = np.argmax(results['like1'])  # best likelihood (objective)
best_params = [results['dec{}'.format(i+1)][best_idx] for i in range(6)]
print("Best parameters:", best_params)

# Evaluate on test dataset
# Prepare test data for simulation
test_setup = BudykoSPOTPY(test)
sim_evap_test, sim_strflow_test = test_setup.simulation(best_params)
obs_evap_test, obs_strflow_test = test_setup.evaluation()

kge_test_str = hydro_kge(sim_strflow_test, obs_strflow_test)
kge_test_ev = hydro_kge(sim_evap_test, obs_evap_test)
print(f"Test KGE Streamflow: {kge_test_str:.4f}")
print(f"Test KGE Evapotranspiration: {kge_test_ev:.4f}")

# Simulate full dataset
full_setup = BudykoSPOTPY(data)
sim_evap_full, sim_strflow_full = full_setup.simulation(best_params)

# Compute final KGE
kge_full_str = hydro_kge(sim_strflow_full, data['hydr_vol'].values)
kge_full_ev = hydro_kge(sim_evap_full, data['Eva_prod_vol'].values)
print(f"Full KGE Streamflow: {kge_full_str:.4f}")
print(f"Full KGE Evapotranspiration: {kge_full_ev:.4f}")
