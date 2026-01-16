import iris,glob
import scipy as sp
from scipy.special import erf
import calendar  # Import the calendar module
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import numpy as np
import iris.coord_systems as cs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.patches as patches
from iris.analysis.cartography import wrap_lons
from matplotlib.colors import LogNorm
import iris.quickplot as qplt
from scipy.spatial import cKDTree
from datetime import datetime, timedelta
import numpy.ma as ma
import pandas as pd
import statistics
import cartopy.feature as cfeature
import iris.plot as iplt
from datetime import datetime, timedelta
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
import metpy.calc as mpcalc
from metpy.units import units
from scipy.stats import norm
import netCDF4 as nc
#from sklearn.metrics import mean_squared_error
from tqdm import tqdm
from netCDF4 import Dataset
import os
from pylab import *


from datetime import datetime, timedelta


#suite = 'do764'
#suite = 'do765'
#suite = 'do766'
suite = 'dn827'

#campaign = 'ATom-1'
#campaign = 'ATom-2'
#campaign = 'ATom-3'
campaign = 'ATom-4'

#num = 1
#num = 2
#num = 3
num = 4

#month = 'aug'
#month = 'feb'
#month = 'oct'
month = 'may'

#year = 2016
#year = 2017
year = 2018


# --- Load all soluble and insoluble mode radii cubes ---
ait_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/aitsol_radius.nc')
acc_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/accsol_radius.nc')
cor_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/corsol_radius.nc')

aitins_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/aitinsol_radius.nc')
accins_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/accinsol_radius.nc')
corins_mode_radius = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/corinsol_radius.nc')

# --- Subset first 11 time steps (or levels) ---
ait_mode_radius = ait_mode_radius[:, 0:11]
acc_mode_radius = acc_mode_radius[:, 0:11]
cor_mode_radius = cor_mode_radius[:, 0:11]

aitins_mode_radius = aitins_mode_radius[:, 0:11]
accins_mode_radius = accins_mode_radius[:, 0:11]
corins_mode_radius = corins_mode_radius[:, 0:11]

# --- Load all soluble and insoluble mode concentration cubes ---
ait_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/aitsol_number.nc')
acc_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/accsol_number.nc')
cor_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/corsol_number.nc')

aitins_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/aitinsol_number.nc')
accins_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/accinsol_number.nc')
corins_mode_conc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/corinsol_number.nc')



# --- Subset first 11 levels (or time steps) ---
ait_mode_conc = ait_mode_conc[:, 0:11]* 1e-6
acc_mode_conc = acc_mode_conc[:, 0:11]* 1e-6
cor_mode_conc = cor_mode_conc[:, 0:11]* 1e-6

aitins_mode_conc = aitins_mode_conc[:, 0:11]* 1e-6
accins_mode_conc = accins_mode_conc[:, 0:11]* 1e-6
corins_mode_conc = corins_mode_conc[:, 0:11]* 1e-6




land_sea_file = '/jet/home/nasch/IMERG_land_sea_mask.nc'  ##just aqua

land_sea_mask = Dataset(land_sea_file, mode='r')
print(land_sea_mask.variables.keys())
land_sea_mask_values = land_sea_mask['landseamask'][:]   ##100 = all water and 0 = all land
land_sea_mask_lat = land_sea_mask['lat'][:]
land_sea_mask_lon = land_sea_mask['lon'][:]


file_path_aerosol = '/jet/home/nasch/tmp_ondemand_ocean_atm200005p_symlink/earhg/observations/ATom_Aerosol_Properties_V2_2111/data/ATom_60s_aerosol_data.nc'  ##just aqua
###All ATom events


aerosol_data = Dataset(file_path_aerosol, mode='r')
print(aerosol_data.variables.keys())
aerosol_data.variables['DateTime_UTC']
#aerosol_data.variables['time']

latitude = aerosol_data.variables['latitude'][:]

# Assuming you have the DateTime_UTC data in a variable
# For demonstration purposes, we'll use a mock variable, replace this with the actual data.
# For example, `DateTime_UTC` data could be extracted like this:
# datetime_utc_data = dataset.variables['DateTime_UTC'][:]
# Here I use a mock dataset for demonstration.

datetime_utc_data = aerosol_data.variables['DateTime_UTC'][:]

# Define the base reference time (1904-01-01 00:00:00)
base_time = datetime(1904, 1, 1)

# Convert seconds to datetime
converted_dates = [base_time + timedelta(seconds=sec) for sec in datetime_utc_data]

aerosol_timestamps = np.array([dt.strftime('%Y-%m-%d %H:%M') for dt in converted_dates])


np.shape(aerosol_timestamps)

latitude_aerosol = aerosol_data.variables['latitude'][:]
longitude_aerosol = aerosol_data.variables['longitude'][:]
num_accum_aerosol = aerosol_data.variables['num_accum'][:]
w_aerosol = aerosol_data.variables['W']
altitude_aerosol = aerosol_data.variables['gps_altitude'][:]




# Define campaign names and corresponding slices
campaign_abbr = ["atom1", "atom2", "atom3", "atom4"]
campaign_slices = [
    slice(0, 5440),       # ATom-1
    slice(5450, 11080),     # ATom-2
    slice(11090, 17810),    # ATom-3
    slice(17820, None)      # ATom-4 (until end)
]

# Dictionary to store variables
data = {}

for campaign_name, campaign_slice in zip(campaign_abbr, campaign_slices):
    # Assign values dynamically with correct slice
    data[f"time_{campaign_name}"] = datetime_utc_data[campaign_slice]
    data[f"latitude_{campaign_name}"] = latitude_aerosol[campaign_slice]
    data[f"longitude_{campaign_name}"] = longitude_aerosol[campaign_slice]
    data[f"num_accum_{campaign_name}"] = num_accum_aerosol[campaign_slice]
    data[f"w_{campaign_name}"] = w_aerosol[campaign_slice]
    data[f"altitude_{campaign_name}"] = altitude_aerosol[campaign_slice]
    

    # Print NaN and masked counts
    print(f"NaN count in num_accum_{campaign_name}:",
          np.isnan(data[f"num_accum_{campaign_name}"]).sum())
    print(f"Masked count in num_accum_{campaign_name}:",
          np.ma.count_masked(data[f"num_accum_{campaign_name}"]))

    # Create a unified mask based on altitude > 1000 or existing mask in num_accum
    mask = (data[f"altitude_{campaign_name}"] > 1000) | np.ma.getmaskarray(data[f"num_accum_{campaign_name}"])
    

    # Apply the mask to all relevant variables
    for key in ["time", "latitude", "longitude", "num_accum", "altitude"]:
        data[f"{key}_{campaign_name}"] = np.ma.masked_where(mask, data[f"{key}_{campaign_name}"])

# If needed, store results as global variables
globals().update(data)




N_frac_ait = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_AIT_1000.nc')
N_frac_ait_ins = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_AITINS_1000.nc')
N_frac_acc = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_ACC_1000.nc')
N_frac_acc_ins = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_ACCINS_1000.nc')
N_frac_cor = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_COR_1000.nc')
N_frac_cor_ins = iris.load_cube(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_CORINS_1000.nc')




time_atom = data[f"time_atom{num}"]
lat_atom = data[f"latitude_atom{num}"]
lon_atom = data[f"longitude_atom{num}"]
alt_atom = data[f"altitude_atom{num}"]
num_accum_atom = data[f"num_accum_atom{num}"]

# --- Convert longitudes from -180:180 to 0:360 ---
lon_atom360 = np.where(lon_atom < 0, lon_atom + 360, lon_atom)

# --- Create a unified mask to filter invalid or out-of-range data ---
mask = (
    np.isnan(lat_atom) |
    np.isnan(lon_atom360) |
    np.isnan(alt_atom) |
    np.isnan(num_accum_atom) ##altitude is already masked >1000 m
)

# --- Apply mask consistently ---
lat_atom_masked = lat_atom[~mask]
lon_atom_masked = lon_atom360[~mask]
alt_atom_masked = alt_atom[~mask]
num_accum_atom_masked = num_accum_atom[~mask]
time_atom_masked = time_atom[~mask]

# --- Convert ATom time (seconds since 1904) to hours since 1970 ---
atom_ref = datetime(1904, 1, 1)
model_ref = datetime(1970, 1, 1)
seconds_offset = (model_ref - atom_ref).total_seconds()
hours_offset = seconds_offset / 3600
time_atom_hours = time_atom_masked / 3600 - hours_offset  # hours since 1970

# --- Model grid coordinates ---
lat_model = acc_mode_conc.coord('latitude').points
lon_model = acc_mode_conc.coord('longitude').points
alt_model = acc_mode_conc.coord('level_height').points
model_time_numeric = acc_mode_conc.coord('time').points  # hours since 1970

# --- Filter ATom points to be within model time range ---
time_in_range_mask = (time_atom_hours >= model_time_numeric.min()) & \
                     (time_atom_hours <= model_time_numeric.max())

lat_atom_masked = lat_atom_masked[time_in_range_mask]
lon_atom_masked = lon_atom_masked[time_in_range_mask]
alt_atom_masked = alt_atom_masked[time_in_range_mask]
num_accum_atom_masked = num_accum_atom_masked[time_in_range_mask]
time_atom_hours = time_atom_hours[time_in_range_mask]

# --- Build KDTree of model lat-lon grid ---
model_grid = np.array([(lat, lon) for lat in lat_model for lon in lon_model])
tree = cKDTree(model_grid)

# --- Prepare ATom observation points ---
obs_points = np.column_stack((lat_atom_masked, lon_atom_masked))

# --- Query nearest model grid point ---
_, indices = tree.query(obs_points)
lat_idx, lon_idx = np.unravel_index(indices, (len(lat_model), len(lon_model)))

# --- Convert masked altitudes to plain array for safe nearest-level search ---
alt_model_arr = alt_model.filled(np.nan)
def find_nearest_level(obs_alt):
    return np.nanargmin(np.abs(alt_model_arr - obs_alt))

level_indices = np.array([find_nearest_level(alt) for alt in alt_atom_masked])

# --- Convert masked model times to plain array for safe nearest-time search ---
model_time_arr = np.array(model_time_numeric)  # in case it's a masked array
def find_nearest_time(obs_time):
    return np.nanargmin(np.abs(model_time_arr - obs_time))

time_indices = np.array([find_nearest_time(t) for t in time_atom_hours])

# --- Convert masked ATom hours back to datetime for inspection ---
atom_times_datetime = np.array([model_ref + timedelta(seconds=h*3600) for h in time_atom_hours])

# --- Print examples ---
print(atom_times_datetime[:20])  # first 20 timestamps
#print(f"Total points after masking: {len(atom_times_datetime)}")
#print(f"Time range: {atom_times_datetime.min()} to {atom_times_datetime.max()}")


print("loading data")
# --- Use core_data to avoid copying large cube data ---
acc_data = acc_mode_conc.data
print("check")
accins_data = accins_mode_conc.data
print("check2")
ait_data = ait_mode_conc.data
aitins_data = aitins_mode_conc.data
cor_data = cor_mode_conc.data
corins_data = corins_mode_conc.data


N_frac_acc_data = N_frac_acc.data
N_frac_acc_ins_data = N_frac_acc_ins.data
N_frac_ait_data = N_frac_ait.data
N_frac_ait_ins_data = N_frac_ait_ins.data
N_frac_cor_data = N_frac_cor.data
N_frac_cor_ins_data = N_frac_cor_ins.data

print("matching data")

# --- Extract only the points corresponding to ATom indices ---
matched_model_accum_values = acc_data[time_indices, level_indices, lat_idx, lon_idx] * \
                              N_frac_acc_data[time_indices, level_indices, lat_idx, lon_idx]
matched_model_accum_insoluble_values = accins_data[time_indices, level_indices, lat_idx, lon_idx] * \
                              N_frac_acc_ins_data[time_indices, level_indices, lat_idx, lon_idx]

matched_model_aitken_values = ait_data[time_indices, level_indices, lat_idx, lon_idx] * \
                              N_frac_ait_data[time_indices, level_indices, lat_idx, lon_idx]

matched_model_aitken_insoluble_values = aitins_data[time_indices, level_indices, lat_idx, lon_idx] * \
                                        N_frac_ait_ins_data[time_indices, level_indices, lat_idx, lon_idx]

matched_model_cor_values = cor_data[time_indices, level_indices, lat_idx, lon_idx] * \
                           N_frac_cor_data[time_indices, level_indices, lat_idx, lon_idx]

matched_model_cor_insoluble_values = corins_data[time_indices, level_indices, lat_idx, lon_idx] * \
                                     N_frac_cor_ins_data[time_indices, level_indices, lat_idx, lon_idx]






# --- Build a dictionary with all columns ---
data_dict = {
    "time": atom_times_datetime,
    "latitude": lat_atom_masked,
    "longitude": lon_atom_masked,
    "altitude_m": alt_atom_masked,
    "num_accum_atom": num_accum_atom_masked,
    "model_accum": matched_model_accum_values,
    "model_accum_insoluble": matched_model_accum_insoluble_values,
    "model_aitken": matched_model_aitken_values,
    "model_aitken_insoluble": matched_model_aitken_insoluble_values,
    "model_cor": matched_model_cor_values,
    "model_cor_insoluble": matched_model_cor_insoluble_values
}

# --- Create a pandas DataFrame ---
df = pd.DataFrame(data_dict)

# --- Save to CSV ---
csv_filename = f"atom_{num}_matched_model.csv"
df.to_csv(csv_filename, index=False)

print(f"CSV saved: {csv_filename}")
print(f"Number of rows: {len(df)}")
print(df.head())