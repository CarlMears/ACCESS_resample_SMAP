# resample_SMAP_to_grid_circular_footprint.py
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from read_SMAP import find_SMAP_L1C_filenames_for_orbit, read_SMAP_L1C_file
from rss_plotting.global_map import plot_global_map
from access_io.access_attr_define import coord_attributes_access 
from access_io.licenses import cc_by_4_license, cc0_license

def output_file_name(*,
                     output_path:Path,
                     input_file_name:str,
                     date_to_do:datetime.date=None):
    
    
    output_path = output_path / date_to_do.strftime('%Y') / date_to_do.strftime('%m') 
    output_path.mkdir(parents=True,exist_ok=True)
    output_file_name = output_path / (input_file_name+'.70km.nc')
    return output_file_name
   
    



def read_smap_wt_file(*,
                      path:Path,
                      target_size:int=70,
                      verbose:bool=True):
    """Reads the L2C weighting files created by resample_SMAP_to_grid_circular_footprint.py"""

    print('Reading L2C Weighting Files:')

    max_locs = 200000
    max_cells_resamp = 1500

    # initialize arrays
    loc_ptr = np.full((720,1440,2),-99,dtype=np.int32)
    weights = np.full((max_locs,max_cells_resamp),np.nan,dtype=np.float64)
    lon_indices = np.full((max_locs,max_cells_resamp),-99,dtype=np.int32)
    lat_indices = np.full((max_locs,max_cells_resamp),-99,dtype=np.int32)
    shapes = np.full((max_locs,2),-99,dtype=np.int32)
    resample_valid = np.zeros((720,1440),dtype=np.int8)
    total_locs = 0

    # step through the target locations
    for itar_lat in range(20,700):
        nc_file = path / f'L2C_weights_{itar_lat:03d}.{target_size:02d}km.v2.nc'
        try:
            ds = xr.open_dataset(nc_file)
        except FileNotFoundError:
            if verbose:
                print(f'File {nc_file} not found')
            continue

        locs = ds['locs'].values
        valid_lats = ds['lat'].values
        valid_lons = ds['lon'].values
        wts = ds['weights'].values
        lat_indices_lat = ds['lat_indices'].values
        lon_indices_lat = ds['lon_indices'].values
        shp = (wts.shape)[2:]
        num_wts = shp[0]*shp[1]
        look = 0
        for loc in locs:
            for look in range(2):
                if np.sum(np.isfinite(wts[loc,look,:,:])) == num_wts:
                    resample_valid[valid_lats[loc],valid_lons[loc]] += 1
                    loc_ptr[valid_lats[loc],valid_lons[loc],look] = total_locs

                    wts_flat = np.reshape(wts[loc,look,:,:],num_wts)
                    weights[total_locs,0:num_wts] = wts_flat

                    lat_indices_flat = np.reshape(lat_indices_lat[loc,:,:],num_wts)
                    lat_indices[total_locs,0:num_wts] = lat_indices_flat

                    lon_indices_flat = np.reshape(lon_indices_lat[loc,:,:],num_wts)
                    lon_indices[total_locs,0:num_wts] = lon_indices_flat
                    
                    shapes[total_locs,:] = shp
                    total_locs += 1

    # truncate arrays to refelct the number of valid locations
    print(f'total_locs = {total_locs}')
    weights = weights[0:total_locs,:]
    lon_indices = lon_indices[0:total_locs,:]
    lat_indices = lat_indices[0:total_locs,:]
    shapes = shapes[0:total_locs,:]

    return (loc_ptr,weights,lon_indices,lat_indices,shapes,resample_valid)

def unflatten_weights(weights, shapes, location):
    shp = shapes[location,:]
    wts = weights[location,0:shp[0]*shp[1]]
    wts = np.reshape(wts,shp)
    return wts

def unwrap_lons(lons):
    """ensure that lons are on the same side of the 0.0/360.0 
    boundary as the center of the grid"""
    shp = lons.shape
    assert len(shp) == 2
    center_loc = [int(shp[0]/2),int(shp[1]/2)]
    center_lon = lons[center_loc[0],center_loc[1]]
    if center_lon > 270.0:
        lons[lons <= 180] = lons[lons <= 180] + 360.0
    
    if center_lon < 90.0:
        lons[lons > 180] = lons[lons > 180] - 360.0
    return lons

if __name__ == '__main__':


    SMAP_weight_path = Path('M:/job_access/resampling/SMAP/weights_for_L2C_files_v3')
    output_path = Path('B:/_access_temp/smap_tb_orbits')
    verbose = False
    make_plots=False

    loc_ptr,weights,lon_indices,lat_indices,shapes,resample_valid = \
            read_smap_wt_file(path=SMAP_weight_path,target_size=70,verbose=verbose)
    print('expanding weights')
    fore_stddev_map = np.full((720,1440),np.nan,dtype=np.float32)   
    aft_stddev_map = np.full((720,1440),np.nan,dtype=np.float32)
    for ilat in range(20,701):
        for ilon in range(0,1440):
            if loc_ptr[ilat,ilon,0] > 0:
                wts = unflatten_weights(weights,shapes,loc_ptr[ilat,ilon,0])
                fore_stddev_map[ilat,ilon] = np.std(wts)
            if loc_ptr[ilat,ilon,1] > 0:
                wts = unflatten_weights(weights,shapes,loc_ptr[ilat,ilon,1])
                aft_stddev_map[ilat,ilon] = np.std(wts)
    
    fig,ax = plot_global_map(fore_stddev_map,
                            title='std. dev. wts, fore',
                            vmin=0,
                            vmax=0.1,
                            cmap='viridis',
                            plt_colorbar=True,
                            plt_coastlines=True) 
    fig,ax = plot_global_map(aft_stddev_map,
                            title='std. dev. wts, aft',
                            vmin=0,
                            vmax=0.1,
                            cmap='viridis',
                            plt_colorbar=True,
                            plt_coastlines=True)
    
    min_std_fore = np.nanmin(fore_stddev_map,axis=1)
    num_bad_fore = np.zeros(720,dtype=np.int32)
    for ilat in range(20,701):
        z = fore_stddev_map[ilat,:]
        z = z[~np.isnan(z)]
        if len(z) > 0:
            num_bad = np.sum(z > 2.0*min_std_fore[ilat])
            num_bad_fore[ilat] = num_bad

    min_std_aft = np.nanmin(aft_stddev_map,axis=1)
    num_bad_aft = np.zeros(720,dtype=np.int32)
    for ilat in range(20,701):
        z = aft_stddev_map[ilat,:]
        z = z[~np.isnan(z)]
        if len(z) > 0:
            num_bad = np.sum(z > 2.0*min_std_aft[ilat])
            num_bad_aft[ilat] = num_bad

    fig3,ax3 = plt.subplots()
    ax3.plot(num_bad_fore)

    fig4,ax4 = plt.subplots()
    ax4.plot(num_bad_aft)
    
    print(num_bad_fore[20:71])
    print(num_bad_aft[20:71])

    print(num_bad_fore[650:701])
    print(num_bad_aft[650:701])

    #plt.show()
    print