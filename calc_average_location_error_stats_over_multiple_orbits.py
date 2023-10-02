# resample_SMAP_to_grid_circular_footprint.py
import numpy as np
import xarray as xr
from pathlib import Path 
import datetime
from read_SMAP import find_SMAP_L1C_filenames_for_orbit,read_SMAP_L1C_file
import matplotlib.pyplot as plt
from rss_gridding.local_earth_grid import LocalEarthGrid
from rss_plotting.global_map import plot_global_map



if __name__ == "__main__":

    verbose = True
    output_path = Path('M:/job_access/resampling/SMAP/weights_for_L2C_files_v2')
    date = datetime.date(2020, 4, 1)
    files = find_SMAP_L1C_filenames_for_orbit(date=date)

    tot_lon_map = np.full((720,1560,2),0.0,dtype=np.float32)
    tot_lat_map = np.full((720,1560,2),0.0,dtype=np.float32)
    tot_sin_eaa_map = np.full((720,1560,2),0.0,dtype=np.float32)
    tot_cos_eaa_map = np.full((720,1560,2),0.0,dtype=np.float32)
    tot_num_map = np.full((720,1560,2),0.0,dtype=np.float32)

    for smap_file in files:
        smap_xr = read_SMAP_L1C_file(path=smap_file)
        print(smap_file)


        #plot some footprint locations:

        smap_lons = smap_xr['cellon'].values
        smap_lats = smap_xr['cellat'].values
        
        smap_lats_ideal = np.round(4.0*(smap_lats+0.125))/4.0 -0.125
        lat_diff = smap_lats_ideal - smap_lats

        smap_lons_ideal = np.round(4.0*(smap_lons+0.125))/4.0 - 0.125
        lon_diff = smap_lons_ideal - smap_lons

        print('lat_diff: ',np.nanmin(lat_diff),np.nanmax(lat_diff))
        print('lon_diff: ',np.nanmin(lon_diff),np.nanmax(lon_diff))
        print('std lat_diff: ',np.nanstd(lat_diff))
        print('std lon_diff: ',np.nanstd(lon_diff))
        print

        lon_diff = np.reshape(lon_diff,(720,1560*2) )
        stddev_lon_by_lat = np.nanstd(lon_diff,axis=1)

        print()
        