import iris,glob
import numpy as np
import iris.plot as iplt








path_m   = '/jet/home/nasch/cylc-run/u-dk832/share/data/History_Data/dk832a.pm2018sep'


lwc_cubes = iris.cube.CubeList()
fraction_cubes = iris.cube.CubeList()
pressure_cubes = iris.cube.CubeList()
theta_cubes = iris.cube.CubeList()




# m01s00i408
list_of_files = [path_m]
print(list_of_files)

for filename in list_of_files:
    lwc = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s00i254'))
    lwc_cubes.append(lwc)
    fraction = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s00i267'))
    fraction_cubes.append(fraction)
    pressure = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s00i408'))
    pressure_cubes.append(pressure)
    theta = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s00i004'))
    theta_cubes.append(theta)





timeseries_lwc = lwc_cubes.concatenate_cube()
timeseries_fraction = fraction_cubes.concatenate_cube()
timeseries_pressure = pressure_cubes.concatenate_cube()
timeseries_theta = theta_cubes.concatenate_cube()






def get_air_density( air_pressure,potential_temperature):
  p0 = iris.coords.AuxCoord(1000.0,
                            long_name='reference_pressure',
                            units='hPa')
  p0.convert_units(air_pressure.units)

  Rd=287.05 # J/kg/K
  cp=1005.46 # J/kg/K
  Rd_cp=Rd/cp

  temperature=potential_temperature*(air_pressure/p0)**(Rd_cp)
  temperature._var_name='temperature'
  R_specific=iris.coords.AuxCoord(287.058,
                                  long_name='R_specific',
                                  units='J-kilogram^-1-kelvin^-1')#J/(kgK)

  air_density=(air_pressure/(temperature*R_specific))
  air_density.long_name='Density of air'
  air_density._var_name='air_density'
  air_density.units='kg m-3'
  temperature.units='K'
  return [air_density, temperature]

[air_density, temperature] = get_air_density(timeseries_pressure, timeseries_theta)







lwc = timeseries_lwc*1000*air_density ##convert from kg/kg to g/m^3
## (kg/kg)*(1000g/kg)*(kg/m^3)  = g/m^3
lwc.units = 'g m-3'
lwc_arr = lwc.data
height = np.array(lwc.coord('level_height').points)

## implement the trapezoid rule to get liquid water path
lwp = np.zeros([144,192])

for j in range(144):
    for k in range(192):
        trap_sum = 0
        for l in range(84):   #85-1
            trap = 0.5*(lwc_arr[l+1,j,k]+lwc_arr[l,j,k])*(height[l+1]-height[l])
            trap_sum += trap

        lwp[j,k] = trap_sum