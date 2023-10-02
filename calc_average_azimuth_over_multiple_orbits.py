# resample_SMAP_to_grid_circular_footprint.py
import numpy as np
import xarray as xr
from pathlib import Path 
import os
import glob
import datetime
from read_SMAP import find_SMAP_L1C_filenames_for_orbit,read_SMAP_L1C_file
from gaussian import Gaussian
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
        smap_eaa  = smap_xr['eaa'].values

        eq_lons_fore = smap_lons[360,:,0]
        longitude_offset = np.round(4.0*(270.0 - eq_lons_fore[1152]))/4.0

        smap_lons = (smap_lons + longitude_offset+720.0)%360.0 
        eq_lons_fore_shifted = smap_lons[360,:,0]
        print(f'longitude_offset = {longitude_offset}')
        fig,ax = plot_global_map(smap_lons[:,:,0], 
                                 title='lons, fore', 
                                 vmin=0, 
                                 vmax=360, 
                                 cmap='twilight',
                                 plt_colorbar=True,
                                 central_longitude=180.0)
        
        fig,ax = plot_global_map(smap_lons[:,:,0]-smap_lons[:,:,1], 
                                 title='lon difference, fore', 
                                 vmin=-.2, 
                                 vmax=.2, 
                                 cmap='BrBG',
                                 plt_colorbar=True,
                                 central_longitude=180.0)
        
        fig2,ax2 = plt.subplots()
        ax2.plot(eq_lons_fore)
        ax2.plot(eq_lons_fore_shifted)

        fig3,ax3 = plt.subplots()
        ax3.plot(smap_lons[690,:,0],marker='+',linestyle='None')
        plt.show()
        sin_smap_eaa = np.sin(np.deg2rad(smap_eaa))
        cos_smap_eaa = np.cos(np.deg2rad(smap_eaa))

        ok = np.all([np.isfinite(smap_lons),
                     np.isfinite(smap_lats),
                     np.isfinite(sin_smap_eaa),
                     np.isfinite(cos_smap_eaa)],
                     axis=0)
        
        tot_lon_map[ok] += smap_lons[ok]
        tot_lat_map[ok] += smap_lats[ok]
        tot_sin_eaa_map[ok] += sin_smap_eaa[ok]
        tot_cos_eaa_map[ok] += cos_smap_eaa[ok]
        tot_num_map[ok] += 1.0

    mean_eaa_map = np.rad2deg(np.arctan2(tot_sin_eaa_map/tot_num_map,tot_cos_eaa_map/tot_num_map))
    
    # fig,ax = plot_global_map(mean_eaa_map[:,:,0], title='eaa, fore', vmin=-180, vmax=180, cmap='twilight',plt_colorbar=True)
    # fig,ax = plot_global_map(mean_eaa_map[:,:,1], title='eaa, aft', vmin=-180, vmax=180, cmap='twilight',plt_colorbar=True)
    # fig,ax = plot_global_map(tot_num_map[:,:,0], title='num, fore', vmin=0.0, vmax=20, cmap='viridis',plt_colorbar=True)
    # fig,ax = plot_global_map(tot_num_map[:,:,1], title='num, aft', vmin=0.0, vmax=20, cmap='viridis',plt_colorbar=True)


    mean_eaa_map_xr = xr.DataArray(mean_eaa_map,dims=['ilat','ilon','look'],coords={'ilat':np.arange(0,720),'ilon':np.arange(0,1560),'look':[0,1]})
    mean_eaa_map_xr.to_netcdf(output_path / 'mean_eaa_map.nc')


    test = xr.open_dataarray(output_path / 'mean_eaa_map.nc')

    print()


        



   