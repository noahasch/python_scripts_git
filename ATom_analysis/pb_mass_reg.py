#!/usr/bin/env python
#qsub -l walltime=00:15:00 -q "shared" -P "clarify" -l select=1:ncpus=1:mem=24GB process_aerosol.py "     

#this script should be run on monsoon xcsc

import numpy as np
import iris,os
print('iris version',iris.__version__)
import iris.fileformats.um
import netCDF4
import copy
import sys
import glob
import datetime
import cf_units as unit
from iris.time import PartialDateTime
import scipy as sp
for arg in sys.argv:
    print(arg)

#job='u-cx608'
#year_month='201802'
#days = ['04']
 
#job='u-cu555'
#year_month='201803'
#days = ['18','20']
file_base=['/ocean/projects/atm200005p/nasch/ATom/simulation_spinup/ATom-4/History_Data/']  ##change if needed
dirs = file_base

print(dirs)
umfile=1
iris_version=2.0


def load_with(filepath,constraint):
  if umfile==1 and iris_version==2:
    from iris.fileformats.um import structured_um_loading
    with structured_um_loading():
      cube = iris.load_cube(filepath,constraint)
  else:
    cube = iris.load_cube(filepath,constraint)
  return cube

def get_air_density( air_pressure,potential_temperature):
  print("Calculating air density...")
  p0 = iris.coords.AuxCoord(1000.0,
                            long_name='reference_pressure',
                            units='hPa')
  p0.convert_units(air_pressure.units)

  Rd=287.05 # J/kg/K
  cp=1005.46 # J/kg/K
  Rd_cp=Rd/cp

  temperature=potential_temperature*(air_pressure/p0)**(Rd_cp)
  #print temperature.data[0,0,0]
  temperature._var_name='temperature'
  R_specific=iris.coords.AuxCoord(287.058,
                                  long_name='R_specific',
                                  units='J-kilogram^-1-kelvin^-1')#J/(kgK)

  air_density=(air_pressure/(temperature*R_specific))
  air_density.long_name='Density of air'
  air_density._var_name='air_density'
  air_density.units='kg m-3'
  temperature.units='K'
  print("Air density calculated.")
  #print air_density.data[0,0,0]
  return [air_density, temperature]


def get_mode_radius(num_cube,mass_list,rho_list,kappa_list, modesigma):
    icube=0
  # Volume mixing ratio in silly units of volume of stuff/kg of air
    mode_vmrs = []
    mode_kappas = []
    for mass_cube in mass_list:
        mode_vmrs.append(mass_cube/rho_list[icube])
        mode_kappas.append(kappa_list[icube]*mode_vmrs[icube])
        icube=icube+1
  # change from volume of stuff per kg of air to per molecule of air
    mode_volume = (28.991e-3/6.022e23)*np.sum(mode_vmrs)
    mode_kappa = (28.991e-3/6.022e23)*np.sum(mode_kappas)/mode_volume
  # volume of stuff per molecule of air divided by num_cube is volume of stuff per particle.
    mode_radius =0.5*iris.analysis.maths.exponentiate((6*mode_volume/num_cube)/(np.pi*np.exp(4.5*(np.log(modesigma))**2.0)),1./3.)
    #print 'radius',mode_radius.data[0,0,0]
    #print 'kappa',mode_kappa.data[0,0,0]
    return [mode_radius, mode_kappa]

def lognormal_cumulative(N,r,rbar,sigma):
    total=(N.data/2)*(1+sp.special.erf(np.log10(r.data/rbar.data)/np.sqrt(2)/np.log10(sigma)))  #number concentration lognormal distribution
    return N.copy(total)                                                                 # below a certain size (dN/dlog*dp)


def process_run_dir(dir):
    print(f"Processing directory: {dir}")
    
    # Include both July and August files
    pp_files = sorted(glob.glob(os.path.join(dir, 'dn827a.pb2018apr')))
    print("Found files:", pp_files)

    for ppfile in pp_files:
        print(f"Processing file: {ppfile}")
        p_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s00i408'))[:,0:11]  # air pressure cube, subset first 11 levels (or time steps)
        theta_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s00i004'))[:,0:11] 
        [air_density, temperature] = get_air_density(p_cube, theta_cube)

        if umfile == 1:
            suffix = ppfile[-7:]
        else:
            suffix = ppfile[-6:-3]

        print("Loading aerosol data...")

        # Load soluble and insoluble number mixing ratios
        ait_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i103'))[:, 0:11]  ##up to 1 km
        acc_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i107'))[:, 0:11]
        cor_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i113'))[:, 0:11]
        aitins_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i119'))[:, 0:11]
        accins_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i122'))[:, 0:11]
        corins_numr_cube = load_with(ppfile, iris.AttributeConstraint(STASH='m01s34i124'))[:, 0:11]
        print("Aerosol data loaded.")

        # Convert to number concentration
        def num_conc(cube):
            return cube * air_density * (6.02e23 / 28.991e-3)

        ait_num_cube = num_conc(ait_numr_cube)
        acc_num_cube = num_conc(acc_numr_cube)
        cor_num_cube = num_conc(cor_numr_cube)
        aitins_num_cube = num_conc(aitins_numr_cube)
        accins_num_cube = num_conc(accins_numr_cube)
        corins_num_cube = num_conc(corins_numr_cube)

        print("Number concentration calculated.")

        # Assign metadata
        cubes = [
            (ait_num_cube, "aitsol_number"),
            (acc_num_cube, "accsol_number"),
            (cor_num_cube, "corsol_number"),
            (aitins_num_cube, "aitinsol_number"),
            (accins_num_cube, "accinsol_number"),
            (corins_num_cube, "corinsol_number"),
        ]

        for cube, name in cubes:
            cube.units = "m-3"
            cube.var_name = name
            cube.long_name = name.replace("_", " ")
            save_path = os.path.join(dir, f"{name}.nc")
            iris.save(cube, save_path)
            print(f"✅ Saved: {save_path}")

        print("All variables saved for this file.\n")



for dir in dirs:
    process_run_dir(dir)

