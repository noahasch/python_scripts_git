## Imports...some won't be needed but I don't feel like sorting through them all

import datetime
import iris,glob
import scipy as sp
from scipy.special import erf
import calendar  # Import the calendar module
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import iris.analysis.cartography as cart
import numpy as np
import iris.coord_systems as cs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from iris.analysis.cartography import wrap_lons
from scipy.spatial import cKDTree
from datetime import datetime, timedelta
import numpy.ma as ma
import pandas as pd
import iris.plot as iplt
import matplotlib.colors as mcolors
import netCDF4 as nc
import os

## load in observational data

lat_naames = lat_met
lon_naames = lon_met
alt_naames = alt_met

##if needed, you can convert obs time to the model format...example below but each case may be different


sec_utc = time_met  # or whatever your column is called
day_start = datetime.datetime(2016, 5, 19, 0, 0, 0, tzinfo=datetime.timezone.utc)
epoch = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
base_hours_1970 = (day_start - epoch).total_seconds() / 3600.0

# convert seconds to hours since 1970
time_naames = base_hours_1970 + sec_utc / 3600.0  ## in model format of hours since 1/1/1970

### If dealing with regional simulations
cube = iris.load_cube('/ocean/projects/atm200005p/nasch/NAAMES/simulations/u-dp878_wvar_0.1_scaled/Regn1/resn_1/RA3/um/aitsol_number_full_run.nc')


##for regional model only!!!!!!!!!!!!
rlat = cube.coord('grid_latitude').points
rlon = cube.coord('grid_longitude').points

rotated_cs = cube.coord_system()

lon_model, lat_model = cart.unrotate_pole(
    rlon, 
    rlat, 
    rotated_cs.grid_north_pole_longitude, 
    rotated_cs.grid_north_pole_latitude
)


### If dealing with global simulations

cube =  iris.load_cube('/jet/home/nasch/cylc-run/u-dp878/share/cycle/20160519T0000Z/glm/um/aitsol_number_full_run_global.nc')

lon_model = cube.coord('longitude').points
lon_model = ((lon_model + 180) % 360) - 180   ##this converts from 0-360 to -180 to 180 which may be needed depending on your obs data
lat_model = cube.coord('latitude').points


###needed for both regional and global simulations

alt_model = cube.coord('level_height').points
model_time_numeric = cube.coord('time').points  


###now find the nearest model indices to observations

model_grid = np.array([(lat, lon) for lat in lat_model for lon in lon_model])
tree = cKDTree(model_grid)

# --- Prepare ATom observation points ---
obs_points = np.column_stack((lat_naames, lon_naames))

# --- Query nearest model grid point ---
_, indices = tree.query(obs_points)
lat_idx, lon_idx = np.unravel_index(indices, (len(lat_model), len(lon_model)))

# --- Convert masked altitudes to plain array for safe nearest-level search ---
# NaN values only seem to be an issue with altitudes from my experience
alt_model_arr = alt_model.filled(np.nan)

def find_nearest_level(obs_alt):
    return np.nanargmin(np.abs(alt_model_arr - obs_alt))

level_indices = np.array([find_nearest_level(alt) for alt in alt_naames])

# --- Convert masked model times to plain array for safe nearest-time search ---
model_time_arr = np.array(model_time_numeric)  # in case it's a masked array

def find_nearest_time(obs_time):
    return np.nanargmin(np.abs(model_time_arr - obs_time))

time_indices = np.array([find_nearest_time(t) for t in time_naames])



###  you can now use lat_idx, lon_idx, level_indices, and time_indices to index into your model data cube

matched_cf = timeseries_cf.data[time_indices, level_indices, lat_idx, lon_idx]