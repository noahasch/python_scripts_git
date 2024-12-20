##note that you need 3 hourly output and n96 resolution for this to work!!!
## 3 hourly mean would be ideal



## add imports
import iris,glob
import numpy as np
import iris.plot as iplt
from tqdm import tqdm
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs

## add file path
path = '/ocean/projects/atm200005p/nasch/CFOG/rose_runs/casim_fixes/casim_on/u-dg148/History_Data/dg148a.pb2018sep'


list_of_files = [path]

print(list_of_files)

## pick the variables you would like to work with...I have a few here
cloud_weights_cubes = iris.cube.CubeList()
cloud_fraction_modis_cubes = iris.cube.CubeList()
optical_depth_modis_cubes = iris.cube.CubeList()
effective_radius_modis_cubes = iris.cube.CubeList()
liquid_water_path_modis_cubes = iris.cube.CubeList()


## these are all COSP model output...not sure how I feel about them but they are meant to replicate satellite data
for filename in list_of_files:
    weights = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s02i330'))
    cloud_weights_cubes.append(weights[:,:,:])
    cloud_fraction = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s02i452'))
    cloud_fraction_modis_cubes.append(cloud_fraction[:,:,:])
    optical_depth = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s02i458'))
    optical_depth_modis_cubes.append(optical_depth[:,:,:])
    radius= iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s02i463'))
    effective_radius_modis_cubes.append(radius[:,:,:])   
    lwp = iris.load_cube(filename,iris.AttributeConstraint(STASH='m01s02i466'))
    liquid_water_path_modis_cubes.append(lwp[:,:,:])



timeseries_weights = cloud_weights_cubes.concatenate_cube()
timeseries_cf = cloud_fraction_modis_cubes.concatenate_cube()
timeseries_tau = optical_depth_modis_cubes.concatenate_cube()
timeseries_eff_rad = effective_radius_modis_cubes.concatenate_cube()
timeseries_lwp = liquid_water_path_modis_cubes.concatenate_cube()

##we need a set of lat and lons to create our new cubes
latitudes = timeseries_weights.coord('latitude').points
longitudes = timeseries_weights.coord('longitude').points



## Aqua satellite matching done here

cf_dict_aqua = {}
weights_dict_aqua = {}
lwp_dict_aqua = {}
tau_dict_aqua = {}
r_eff_dict_aqua = {}


## this next bit is tricky


## we break the earth down into 8 equal parts
## and for aqua we want to sync up the 1/8th part of the domain
## where the local time is ~1:30 PM (13:30)
## I do this by finding the range of 24 longitude points out of the 192 
## for each 3 hour time interval that match the satellite
## note that model longitude increasing in the opposite direction to the satellite rotation
## so the values of longitude get smaller as we move west from the 
## antimerdian towards the prime meridian, following the satellite path
## as the model index is 0 at the prime meridian, 96 at anti (moving east), and 192 back at prime meridian

## tqdm just shows how long the for loop has to run, which can take a while
for i in tqdm(range(0,8)):
    if i < 4: 
        cf_data = timeseries_cf[i::8, :, 108-((i+1)*24):108-(i*24)].data 
        weights_data = timeseries_weights[i::8, :, 108-((i+1)*24):108-(i*24)].data
        lwp_data = timeseries_lwp[i::8, :, 108-((i+1)*24):108-(i*24)].data

        tau_data = timeseries_tau[i::8, :, 108-((i+1)*24):108-(i*24)].data
        r_eff_data = timeseries_eff_rad[i::8, :, 108-((i+1)*24):108-(i*24)].data
        lons = iris.coords.DimCoord(longitudes[108-((i+1)*24):108-(i*24)], standard_name='longitude', units='degrees_east')

    elif i == 4:  # For handling data across the prime meridian
        cf_data = np.concatenate((timeseries_cf[i::8, :, -12:].data, timeseries_cf[i::8, :, :12].data), axis = 2)
        weights_data = np.concatenate((timeseries_weights[i::8, :, -12:].data, timeseries_weights[i::8, :, :12].data), axis = 2)
        lwp_data = np.concatenate((timeseries_lwp[i::8, :, -12:].data, timeseries_lwp[i::8, :, :12].data), axis = 2)

        tau_data = np.concatenate((timeseries_tau[i::8, :, -12:].data, timeseries_tau[i::8, :, :12].data), axis = 2)
        r_eff_data = np.concatenate((timeseries_eff_rad[i::8, :, -12:].data, timeseries_eff_rad[i::8, :, :12].data), axis = 2)
        
        # Concatenate and shift longitudes to ensure they are monotonic
        lon_concat = np.concatenate((longitudes[-12:], longitudes[:12]))
        lon_concat[lon_concat > 180] -= 360  # Shift longitudes >180 to negative values
        
        lons = iris.coords.DimCoord(lon_concat, standard_name='longitude', units='degrees_east')

    else:
        cf_data = timeseries_cf[i::8,:,180-((i-4)*24):180-((i-5)*24)].data  ##do i+4 for aqua and i+3 for terra
        weights_data = timeseries_weights[i::8,:,180-((i-4)*24):180-((i-5)*24)].data
        lwp_data = timeseries_lwp[i::8,:,180-((i-4)*24):180-((i-5)*24)].data

        tau_data = timeseries_tau[i::8,:,180-((i-4)*24):180-((i-5)*24)].data
        r_eff_data = timeseries_eff_rad[i::8,:,180-((i-4)*24):180-((i-5)*24)].data

        lons = iris.coords.DimCoord(longitudes[180-((i-4)*24):180-((i-5)*24)], standard_name='longitude', units='degrees_east')  

    time = iris.coords.DimCoord(timeseries_weights.coord('time').points[i::8], standard_name='time', units='hours since 1970-01-01 00:00:00')
    lats = iris.coords.DimCoord(latitudes, standard_name='latitude', units='degrees_north')
    
    # Create cubes
    cf = iris.cube.Cube(cf_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    weights = iris.cube.Cube(weights_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    lwp = iris.cube.Cube(lwp_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    tau = iris.cube.Cube(tau_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    r_eff = iris.cube.Cube(r_eff_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    
    # Store in dictionaries
    cf_dict_aqua[f'cf_{i+1}'] = cf
    weights_dict_aqua[f'weights_{i+1}'] = weights
    lwp_dict_aqua[f'lwp_{i+1}'] = lwp
    tau_dict_aqua[f'tau_{i+1}'] = tau
    r_eff_dict_aqua[f'r_eff_{i+1}'] = r_eff


## now we need to convert our dictionaries (with a temporal dimension) into a temporal mean
## pick the amount of days you want to use!!
## my simulation was for three months so I can select of range to isolate single months


days = np.arange(0, 30)  #September
#days = np.arange(30, 61)  #October
#days = np.arange(61,91)  #November


## intialize arrays
lwp_in_cloud_step = []
lwp_in_cloud_aqua = np.zeros((144,192))
lwp_grid_av_aqua = np.zeros((144,192))
nd_in_cloud_aqua = np.zeros((144,192))
cf_list_aqua = np.zeros((144,192))

sqrt_5 = np.sqrt(5)
two_pi_08 = 2 * np.pi * 0.8
constant_factor = (sqrt_5 / two_pi_08) * ((0.7 * 1.4 * 10**-6) / (2 * 1000))**0.5
for i in range(0,8):
    tau_data = tau_dict_aqua[f'tau_{i+1}'][days[0]:days[-1]+1].data  ##day-1 to offset first day being 0th index position
    cf_data = cf_dict_aqua[f'cf_{i+1}'][days[0]:days[-1]+1].data
    rad_data = r_eff_dict_aqua[f'r_eff_{i+1}'][days[0]:days[-1]+1].data
    lwp_data = lwp_dict_aqua[f'lwp_{i+1}'][days[0]:days[-1]+1].data*10**3  ##kg/m^2 to g/m^2
    

    tau_mask = ma.masked_where(tau_data < 0, tau_data)        #mask first and then take the average
    cf_mask = ma.masked_where(cf_data < 0.1, cf_data)        # masking meant to replicate grosvenor 2018 but really closer to Quaas 2006 
    rad_mask = ma.masked_where(rad_data < 0 * 10**-6, rad_data) #refer to this paper: https://amt.copernicus.org/articles/15/3875/2022/
    lwp_mask = ma.masked_where(lwp_data < 10, lwp_data)
  

    tau_in_cloud = tau_data / cf_mask   #masking before I calculate in cloud
    rad_in_cloud = rad_data / cf_mask   #may need to think more about this
    lwp_in_c  = lwp_mask / cf_mask


    # Use broadcasting to perform the calculation on the entire array
    Nd_list = constant_factor * (tau_in_cloud / rad_in_cloud**5)**0.5 * 10**-6 ##using Nd equation in https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2017RG000593 
    #grid averaged values here
    mean_lwp = np.nanmean(lwp_data, axis = 0)  
    mean_cf = np.nanmean(cf_dict_aqua[f'cf_{i+1}'][0:len(days)].data, axis = 0)
    mean_nd = np.nanmean(Nd_list , axis = 0)
    
    #in cloud
    mean_lwp_in_cloud = np.nanmean(lwp_in_c, axis = 0)  

    ##assign the values from the dictionaries to the correct indexes in the zeros arrays
    ##because i=0 corresponds with the antimerdian, we must put these 24 points centered about index 96
    if i < 4:
        lwp_in_cloud_aqua[:,108-((i+1)*24):108-(i*24)] = mean_lwp_in_cloud[:,0:24]  #put this back in the correct location so it matches with lon list below
        lwp_grid_av_aqua[:,108-((i+1)*24):108-(i*24)] = mean_lwp[:,0:24]
        nd_in_cloud_aqua[:,108-((i+1)*24):108-(i*24)] = mean_nd[:,0:24]
        cf_list_aqua[:,108-((i+1)*24):108-(i*24)] = mean_cf[:,0:24]
        
    ##need to split dictionary values at the prime merdian as model longitude array starts and ends here    
    elif i == 4:
        #print(np.shape(mean_lwp))
        lwp_in_cloud_aqua[:,180:192] = mean_lwp_in_cloud[:,0:12]  #put this back in the correct location so it matches with lon list below
        lwp_in_cloud_aqua[:,0:12] = mean_lwp_in_cloud[:,12:24]

        lwp_grid_av_aqua[:,180:192] = mean_lwp[:,0:12] #put this back in the correct location so it matches with lon list below
        lwp_grid_av_aqua[:,0:12] = mean_lwp[:,12:24]

        nd_in_cloud_aqua[:,180:192] = mean_nd[:,0:12] #put this back in the correct location so it matches with lon list below
        nd_in_cloud_aqua[:,0:12] = mean_nd[:,12:24]

        cf_list_aqua[:,180:192] = mean_cf[:,0:12]
        cf_list_aqua[:,0:12] = mean_cf[:,12:24]



    else:
        lwp_in_cloud_aqua[:,180-((i-4)*24):180-((i-5)*24)] = mean_lwp_in_cloud[:,0:24]
        lwp_grid_av_aqua[:,180-((i-4)*24):180-((i-5)*24)] = mean_lwp[:,0:24]
        nd_in_cloud_aqua[:,180-((i-4)*24):180-((i-5)*24)] = mean_nd[:,0:24]  #put this back in the correct location so it matches with lon list below  
        cf_list_aqua[:,180-((i-4)*24):180-((i-5)*24)] = mean_cf[:,0:24]







#### now we must do the same thing for terra


cf_dict_terra = {}
weights_dict_terra = {}
lwp_dict_terra = {}
tau_dict_terra = {}
r_eff_dict_terra = {}


### Note that everything is shifted to the left by one model timestep
## This is because terra passes at 10:30 local time, and thus reaches prime meridian one timestep earlier
for i in tqdm(range(0,8)):
    if i < 3:  # First half
        cf_data = timeseries_cf[i::8, :, 84-((i+1)*24):84-(i*24)].data  # do i+4 for Aqua and i+3 for Terra
        weights_data = timeseries_weights[i::8, :, 84-((i+1)*24):84-(i*24)].data
        lwp_data = timeseries_lwp[i::8, :, 84-((i+1)*24):84-(i*24)].data

        tau_data = timeseries_tau[i::8, :, 84-((i+1)*24):84-(i*24)].data
        r_eff_data = timeseries_eff_rad[i::8, :, 84-((i+1)*24):84-(i*24)].data
        lons = iris.coords.DimCoord(longitudes[84-((i+1)*24):84-(i*24)], standard_name='longitude', units='degrees_east')

    elif i == 3:  # For handling data across the prime meridian
        cf_data = np.concatenate((timeseries_cf[i::8, :, -12:].data, timeseries_cf[i::8, :, :12].data), axis = 2)
        weights_data = np.concatenate((timeseries_weights[i::8, :, -12:].data, timeseries_weights[i::8, :, :12].data), axis = 2)
        lwp_data = np.concatenate((timeseries_lwp[i::8, :, -12:].data, timeseries_lwp[i::8, :, :12].data), axis = 2)

        tau_data = np.concatenate((timeseries_tau[i::8, :, -12:].data, timeseries_tau[i::8, :, :12].data), axis = 2)
        r_eff_data = np.concatenate((timeseries_eff_rad[i::8, :, -12:].data, timeseries_eff_rad[i::8, :, :12].data), axis = 2)
        
        # Concatenate and shift longitudes to ensure they are monotonic
        lon_concat = np.concatenate((longitudes[-12:], longitudes[:12]))
        lon_concat[lon_concat > 180] -= 360  # Shift longitudes >180 to negative values
        
        lons = iris.coords.DimCoord(lon_concat, standard_name='longitude', units='degrees_east')

    else:
        cf_data = timeseries_cf[i::8,:,180-((i-3)*24):180-((i-4)*24)].data  ##do i+4 for aqua and i+3 for terra
        weights_data = timeseries_weights[i::8,:,180-((i-3)*24):180-((i-4)*24)].data
        lwp_data = timeseries_lwp[i::8,:,180-((i-3)*24):180-((i-4)*24)].data

        tau_data = timeseries_tau[i::8,:,180-((i-3)*24):180-((i-4)*24)].data
        r_eff_data = timeseries_eff_rad[i::8,:,180-((i-3)*24):180-((i-4)*24)].data

        lons = iris.coords.DimCoord(longitudes[180-((i-3)*24):180-((i-4)*24)], standard_name='longitude', units='degrees_east')  

    time = iris.coords.DimCoord(timeseries_weights.coord('time').points[i::8], standard_name='time', units='hours since 1970-01-01 00:00:00')
    lats = iris.coords.DimCoord(latitudes, standard_name='latitude', units='degrees_north')
    
    # Create cubes
    cf = iris.cube.Cube(cf_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    weights = iris.cube.Cube(weights_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    lwp = iris.cube.Cube(lwp_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    tau = iris.cube.Cube(tau_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    r_eff = iris.cube.Cube(r_eff_data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(time, 0), (lats, 1), (lons, 2)])
    
    # Store in dictionaries
    cf_dict_terra[f'cf_{i+1}'] = cf
    weights_dict_terra[f'weights_{i+1}'] = weights
    lwp_dict_terra[f'lwp_{i+1}'] = lwp
    tau_dict_terra[f'tau_{i+1}'] = tau
    r_eff_dict_terra[f'r_eff_{i+1}'] = r_eff



days = np.arange(0, 30)  #September
#days = np.arange(30, 61)  #October
#days = np.arange(61,91)  #November

lwp_in_cloud_step = []
lwp_in_cloud_terra = np.zeros((144,192))
lwp_grid_av_terra = np.zeros((144,192))
nd_in_cloud_terra = np.zeros((144,192))

cf_list_terra = np.zeros((144,192))

sqrt_5 = np.sqrt(5)
two_pi_08 = 2 * np.pi * 0.8
constant_factor = (sqrt_5 / two_pi_08) * ((0.7 * 1.4 * 10**-6) / (2 * 1000))**0.5
for i in range(0,8):
    tau_data = tau_dict_terra[f'tau_{i+1}'][days[0]:days[-1]+1].data  ##day-1 to offset first day being 0th index position
    cf_data = cf_dict_terra[f'cf_{i+1}'][days[0]:days[-1]+1].data
    rad_data = r_eff_dict_terra[f'r_eff_{i+1}'][days[0]:days[-1]+1].data
    lwp_data = lwp_dict_terra[f'lwp_{i+1}'][days[0]:days[-1]+1].data*10**3
    

    tau_mask = ma.masked_where(tau_data < 0, tau_data)        #mask first and then take the average
    cf_mask = ma.masked_where(cf_data <= 0.1, cf_data)        # masking meant to replicate grosvenor 2018 bu really closer to Quaas 2006 
    rad_mask = ma.masked_where(rad_data < 0 * 10**-6, rad_data)
    lwp_mask = ma.masked_where(lwp_data < 10, lwp_data)

    tau_in_cloud = tau_data / cf_mask
    rad_in_cloud = rad_data / cf_mask
    lwp_in_c  = lwp_mask / cf_mask


    # Use broadcasting to perform the calculation on the entire array
    Nd_list = constant_factor * (tau_in_cloud / rad_in_cloud**5)**0.5 * 10**-6  

    #grid averaged
    mean_lwp = np.nanmean(lwp_data, axis = 0)  
    mean_cf = np.nanmean(cf_dict_terra[f'cf_{i+1}'][days[0]:days[-1]+1].data, axis = 0)
    mean_nd = np.nanmean(Nd_list , axis = 0)

    ##in cloud
    mean_lwp_in_cloud = np.nanmean(lwp_in_c, axis = 0)  #in cloud


    if i < 3:
        lwp_in_cloud_terra[:,84-((i+1)*24):84-(i*24)] = mean_lwp_in_cloud[:,0:24]  #put this back in the correct location so it matches with lon list below
        lwp_grid_av_terra[:,84-((i+1)*24):84-(i*24)] = mean_lwp[:,0:24]
        nd_in_cloud_terra[:,84-((i+1)*24):84-(i*24)] = mean_nd[:,0:24]
        cf_list_terra[:,84-((i+1)*24):84-(i*24)] = mean_cf[:,0:24]
        
        
    elif i == 3:
        #print(np.shape(mean_lwp))
        lwp_in_cloud_terra[:,180:192] = mean_lwp_in_cloud[:,0:12]  #put this back in the correct location so it matches with lon list below
        lwp_in_cloud_terra[:,0:12] = mean_lwp_in_cloud[:,12:24]

        lwp_grid_av_terra[:,180:192] = mean_lwp[:,0:12] #put this back in the correct location so it matches with lon list below
        lwp_grid_av_terra[:,0:12] = mean_lwp[:,12:24]

        nd_in_cloud_terra[:,180:192] = mean_nd[:,0:12] #put this back in the correct location so it matches with lon list below
        nd_in_cloud_terra[:,0:12] = mean_nd[:,12:24]

        cf_list_terra[:,180:192] = mean_cf[:,0:12]
        cf_list_terra[:,0:12] = mean_cf[:,12:24]



    else:
        lwp_in_cloud_terra[:,180-((i-3)*24):180-((i-4)*24)] = mean_lwp_in_cloud[:,0:24]
        lwp_grid_av_terra[:,180-((i-3)*24):180-((i-4)*24)] = mean_lwp[:,0:24]
        nd_in_cloud_terra[:,180-((i-3)*24):180-((i-4)*24)] = mean_nd[:,0:24]  #put this back in the correct location so it matches with lon list below  
        cf_list_terra[:,180-((i-3)*24):180-((i-4)*24)] = mean_cf[:,0:24]


##### for model output you can use the following code 

figures = []
plt.ioff()



data = (lwp_in_cloud_aqua[20:-20]+lwp_in_cloud_terra[20:-20])/2
mean = np.nanmean(data)
print(mean)
data = ma.masked_where(data == 0, data) 

## I cut off the poles as satellites have issues with sharp zenith angles

lats = iris.coords.DimCoord(latitudes[20:-20], standard_name='latitude', units='degrees_north')
lons = iris.coords.DimCoord(longitudes, standard_name='longitude', units='degrees_east')
cube = iris.cube.Cube(data, var_name='my_variable', units='cm-3', dim_coords_and_dims=[(lats, 0), (lons, 1)])
cube.rename('Eff_rad')

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())

out = iplt.pcolormesh(cube, cmap='jet', vmin = 0, vmax = 500)


ax.add_feature(cfeature.OCEAN, facecolor='#CCFEFF')
ax.add_feature(cfeature.LAKES, facecolor='#CCFEFF')
ax.add_feature(cfeature.RIVERS, facecolor='#CCFEFF')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.LAND, facecolor='#FFE9B5')

# Add colorbar
cbar = plt.colorbar(out, ax=ax, orientation='vertical', pad=0.05, aspect=10, shrink = .5)
cbar.set_label('g/m$^2$', fontsize = 14)
cbar.ax.tick_params(labelsize = 12)

plt.title(f'COSP LWP In Cloud PC2 Sep 2018 (A+T) mean = {mean:.1f} g/m$^2$ ')

# Save the figure to a file
plt.savefig('cosp_lwp_in_cloud_pc2_sep2018.png', dpi=300, bbox_inches='tight')
