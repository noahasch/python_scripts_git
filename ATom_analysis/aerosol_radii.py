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

#do764a.pb2016jul
#do765a.pb2017jan

num = 4
path_b = f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/ATom-{num}/History_Data/dn827a.pb2018apr'
#path_b = '/jet/home/nasch/cylc-run/u-dn188/share/data/History_Data/dn188a.pb2018jan'

list_of_files = [path_b]

print(list_of_files)

ait_sol_cubes = iris.cube.CubeList()
acc_sol_cubes = iris.cube.CubeList()
cor_sol_cubes = iris.cube.CubeList()
ait_insol_cubes = iris.cube.CubeList()
acc_insol_cubes = iris.cube.CubeList()
cor_insol_cubes = iris.cube.CubeList()


for filename in list_of_files:
    ait_sol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i402'))
    ait_sol_cubes.append(ait_sol)
    acc_sol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i403'))
    acc_sol_cubes.append(acc_sol)
    cor_sol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i404'))
    cor_sol_cubes.append(cor_sol)

    ait_insol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i405'))
    ait_insol_cubes.append(ait_insol)
    acc_insol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i406'))
    acc_insol_cubes.append(acc_insol)
    cor_insol = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s38i407'))
    cor_insol_cubes.append(cor_insol)


aitsol_diameter = ait_sol_cubes.concatenate_cube()
aitsol_radius = aitsol_diameter / 2.0
accsol_diameter = acc_sol_cubes.concatenate_cube()
accsol_radius = accsol_diameter / 2.0
corsol_diameter = cor_sol_cubes.concatenate_cube()
corsol_radius = corsol_diameter / 2.0

aitinsol_diameter = ait_insol_cubes.concatenate_cube()
aitinsol_radius = aitinsol_diameter / 2.0
accinsol_diameter = acc_insol_cubes.concatenate_cube()
accinsol_radius = accinsol_diameter / 2.0
corinsol_diameter = cor_insol_cubes.concatenate_cube()
corinsol_radius = corinsol_diameter / 2.0

# --- Define base directory ---
base_dir = f'/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/ATom-{num}/History_Data/'

# --- Assign metadata and save each cube ---
for cube, name in [
    (aitsol_radius, 'aitsol_radius'),
    (accsol_radius, 'accsol_radius'),
    (corsol_radius, 'corsol_radius'),
    (aitinsol_radius, 'aitinsol_radius'),
    (accinsol_radius, 'accinsol_radius'),
    (corinsol_radius, 'corinsol_radius'),
]:
    cube.units = 'm'
    cube.var_name = name
    cube.long_name = name.replace('_', ' ')
    
    # Save each cube individually
    save_path = os.path.join(base_dir, f'{name}.nc')
    iris.save(cube, save_path)
    print(f"✅ Saved: {save_path}")

print("All concatenated mode radius cubes saved successfully!")
