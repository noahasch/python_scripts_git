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


#suite = 'do764'
#suite = 'do765'
#suite = 'do766'
suite = 'dn827'

#campaign = 'ATom-1'
#campaign = 'ATom-2'
#campaign = 'ATom-3'
campaign = 'ATom-4'

#month = 'jul'
#month = 'jan'
#month = 'sep'
month = 'apr'

#year = 2016
#year = 2017
year = 2018


EPS = 1e-20
file_base=[f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/']
dirs = file_base

print(dirs)
umfile=1
iris_version=2.0

#def lognormal_cumulative(N,r,rbar,sigma):
#    total=(N/2)*(1+sp.special.erf(np.log10(r/rbar)/np.sqrt(2)/np.log10(sigma)))  #number concentration lognormal distribution
#    return total

def load_with(filepath,constraint):
  if umfile==1 and iris_version==2:
    from iris.fileformats.um import structured_um_loading
    with structured_um_loading():
      cube = iris.load_cube(filepath,constraint)
  else:
    cube = iris.load_cube(filepath,constraint)
  return cube     





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


pp_files=sorted(glob.glob(f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/{suite}a.pm{year}{month}') )

print(pp_files)
ait_numr_cube = load_with(pp_files,iris.AttributeConstraint(STASH='m01s34i103'))
latitudes = ait_numr_cube.coord('latitude').points
longitudes = ait_numr_cube.coord('longitude').points
height = ait_numr_cube.coord('level_height').points[:11]


def lognormal_cumulative_cube(N_cube, r_threshold, r_cube, sigmas):
    """
    Apply lognormal cumulative to all elements of iris cubes.
    
    Parameters:
    - N_cube: Iris cube of number concentrations
    - r_threshold: Scalar or cube of threshold radius
    - r_cube: Iris cube of mode radius
    - sigma: Geometric standard deviation
    - above: If True, cumulative for r > threshold; if False, r < threshold
    
    Returns:
    - cube_out: Iris cube with cumulative number above/below threshold
    """
    data_mode = N_cube.core_data()
    rad_data = r_cube.core_data()
    erf_arg = np.log10(r_threshold / rad_data) / (np.sqrt(2) * np.log10(sigma))

    data = (data_mode / 2.0) * (1 + sp.special.erf(erf_arg))

    return N_cube.copy(data)


# --- Define cubes in size order: AIT, ACC, COR ---
number_conc_list = [
    [ait_mode_conc, aitins_mode_conc],
    [acc_mode_conc, accins_mode_conc],
    [cor_mode_conc, corins_mode_conc]
]

radius_list = [
    [ait_mode_radius, aitins_mode_radius],
    [acc_mode_radius, accins_mode_radius],
    [cor_mode_radius, corins_mode_radius]
]

sigma_values_list = [
    [1.59, 1.59],  # Aitken
    [1.40, 1.59],  # Accumulation
    [2.0, 1.59]    # Coarse
]

name_list_all = [
    ['AIT', 'AITINS'],
    ['ACC', 'ACCINS'],
    ['COR', 'CORINS']
]

# --- Loop through sets ---
for set_idx in range(len(number_conc_list)):
    number_conc = number_conc_list[set_idx]
    radius = radius_list[set_idx]
    sigma_values = sigma_values_list[set_idx]
    name_list = name_list_all[set_idx]

    for i in range(len(number_conc)):
        num_cube = number_conc[i]
        rad_cube = radius[i]
        sigma = sigma_values[i]
        name = name_list[i]

        if 'AIT' in name:
            # r > 40 nm
            N_below_40 = lognormal_cumulative_cube(num_cube, 40e-9, rad_cube, sigma)
            frac_arr_cube = (num_cube.core_data() - N_below_40.core_data()) / (num_cube.core_data())

        elif 'ACC' in name:
            # Double cutoff: 40 nm < r < 500 nm
            N_below_40 = lognormal_cumulative_cube(num_cube, 40e-9, rad_cube, sigma)
            N_below_500 = lognormal_cumulative_cube(num_cube, 500e-9, rad_cube, sigma)
            frac_arr_cube = (N_below_500.core_data() - N_below_40.core_data()) / (num_cube.core_data())

        elif 'COR' in name:
            # r < 500 nm
            N_below_500 = lognormal_cumulative_cube(num_cube, 500e-9, rad_cube, sigma)
            frac_arr_cube = N_below_500.core_data() / (num_cube.core_data())

        # Create new cube
        cube_frac = num_cube.copy(frac_arr_cube)
        cube_frac.long_name = f"Fraction for {name}"
        cube_frac.var_name = f"N_frac_{name}"

        print(f"Cube for {name}")
        print(cube_frac)
        print("Mean fraction:", np.mean(cube_frac.core_data()))

        # Save the cube
        iris.save(cube_frac, f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/{campaign}/History_Data/N_frac_{name}_1000.nc')


        print(f"Saved N_frac_{name}_1000.nc")
