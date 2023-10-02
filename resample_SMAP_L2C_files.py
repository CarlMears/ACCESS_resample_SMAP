# resample_SMAP_to_grid_circular_footprint.py
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from read_SMAP import find_SMAP_L1C_filenames_for_orbit, read_SMAP_L1C_file
from rss_plotting.global_map import plot_global_map
from access_io.access_attr_define import coord_attributes_access 
from access_io.cc_by_4_license_string import cc_by_4_license_string

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
    import argparse

    parser = argparse.ArgumentParser(   
        description=("Resample SMAP L2C files to 70 km grid")
    )

    parser.add_argument(
        "--smap_wt_root", type=Path, help="Root directory SMAP resampled weights"
    )
    parser.add_argument(
        "--output_root", type=Path, help="Root directory store output files"
    )
    parser.add_argument(
        "--start_date",
        type=datetime.date.fromisoformat,
        help="First Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "--end_date",
        type=datetime.date.fromisoformat,
        help="Last Day to process, as YYYY-MM-DD",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print verbose output"
    )
    parser.add_argument(
        "--make_plots", action="store_true", help="Make plots"
    )

    args = parser.parse_args()

    SMAP_weight_path: Path = args.smap_wt_root
    output_path: Path = args.output_root

    start_date = args.start_date
    end_date = args.end_date 

    verbose = args.verbose
    make_plots = args.make_plots

    # SMAP_weight_path = Path('M:/job_access/resampling/SMAP/weights_for_L2C_files_v2')
    # output_path = Path('B:/_access_temp/smap_tb_orbits')
    # verbose = False
    # make_plots=False

    loc_ptr,weights,lon_indices,lat_indices,shapes,resample_valid = \
            read_smap_wt_file(path=SMAP_weight_path,target_size=70,verbose=verbose)

    if make_plots:
        fig,ax = plot_global_map(resample_valid, title='resample_valid', vmin=0, vmax=2, cmap='viridis')  # noqa: E501



    date_to_do = start_date

    while date_to_do <= end_date:  

        files = find_SMAP_L1C_filenames_for_orbit(date=date_to_do)

        for file_to_process in files:

            output_file_nc = output_file_name(output_path=output_path,
                                            input_file_name=file_to_process.stem,
                                            date_to_do=date_to_do)

            smap_xr = read_SMAP_L1C_file(path=file_to_process)

            list_of_keys_to_copy = ['Conventions',
                                    'processing_level',
                                    'level',
                                    'resolution',
                                    'history',
                                    'date_created',
                                    'institution',
                                    'source',
                                    'platform',
                                    'instrument',
                                    'keywords',
                                    'keywords_vocabulary',
                                    'standard_name_vocabulary',
                                    'creator_name',
                                    'creator_email',
                                    'creator_url',
                                    'publisher_name',
                                    'publisher_email',
                                    'publisher_url',
                                    'id',
                                    'naming_authority',
                                    'dataset_citation_authors',
                                    'dataset_citation_year',
                                    'dataset_citation_product',
                                    'dataset_citation_version',
                                    'dataset_citation_institution',
                                    'dataset_citation_url',
                                    'netCDF_version_id',
                                    'comment',
                                    'references',
                                    'time_coverage_start',
                                    'time_coverage_end',
                                    'orbit_number',
                                    'time_coverage_duration',
                                    'time_coverage_resolution',
                                    'cdm_data_type',
                                    'geospatial_bounds',
                                    'geospatial_lat_min',
                                    'geospatial_lat_max',
                                    'geospatial_lat_resolution',
                                    'geospatial_lat_units',
                                    'geospatial_lon_min',
                                    'geospatial_lon_max',
                                    'geospatial_lon_resolution',
                                    'geospatial_lon_units',
                                    'geospatial_bounds_vertical_crs',
                                    'geospatial_vertical_min',
                                    'geospatial_vertical_max'
                                    ]
            output_attrs = {key:smap_xr.attrs[key] for key in list_of_keys_to_copy}
            output_attrs['title'] = 'SMAP L2C Radiance resampled to 70 km grid'
            output_attrs['version'] = 'V1.0'
            output_attrs['summary'] = 'SMAP L2C Radiance resampled to 70 km grid'
            output_attrs['date_created'] = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
            output_attrs['project'] = 'NASA ACCESS-019'
            output_attrs['keywords'] = 'SMAP, RADIANCE, RESAMPLED, NASA, JPL, RSS'
            output_attrs['license'] = cc_by_4_license_string()
            
            lons = smap_xr['cellon'].values
            lats = smap_xr['cellat'].values
            tb = smap_xr['tb_toa'].values
            time = smap_xr['time'].values
            azim = smap_xr['eaa'].values
            inc = smap_xr['eia'].values

            output_map = np.full((721,1440,2,4),np.nan,dtype=np.float32)   
            time_map = np.full((721,1440,2),np.nan,dtype=np.float64)   
            azim_map = np.full((721,1440,2),np.nan,dtype=np.float32)
            inc_map = np.full((721,1440,2),np.nan,dtype=np.float32)

            stddev_wts = np.full((721,1440,2),np.nan,dtype=np.float32)

            min_ilon = 400
            print(f'Processing {file_to_process}')
            for ilat_out in range(20,700):
                #print(f'ilat_out = {ilat_out}')
                if ilat_out%10 == 0:
                    print('.',end='')
                for ilon_out in range(0,1440):
                   
                    for look in range(2):
                        #find the location in the list of wts, etc for this output lat/lon location
                        loc_index = loc_ptr[ilat_out,ilon_out,look]

                        if loc_index >= 0: #not valid is -99
                            # if ilon_out < 2:
                            #     print(f'loc_index = {loc_index}')
                            #     print
                            wts = unflatten_weights(weights, shapes, loc_index)
                            ilats = unflatten_weights(lat_indices, shapes, loc_index)
                            ilons = unflatten_weights(lon_indices, shapes, loc_index)
                            stddev_wts[ilat_out,ilon_out,look] = np.std(wts)
                            if np.std(wts) > 0.1:
                                plt.imshow(wts)
                                print(f'std dev wts = {np.std(wts)}')
                                print(f'ilon_out = {ilon_out}, ilat_out = {ilat_out}')
                                continue
                            
                            # lons_subset = input_lon_grid[ilats,ilons]
                            lons_subset2 = lons[:,:,look][ilats,ilons]
                            lons_subset2 = unwrap_lons(lons_subset2)
                            
                            #lats_subset = input_lat_grid[ilats,ilons]
                            lats_subset2 = lats[:,:,look][ilats,ilons]
                            time_subset2 = time[:,:,look][ilats,ilons]
                            azim_subset2 = azim[:,:,look][ilats,ilons]
                            # if ilats[5,5] > 200:
                            #     if (np.nanmax(azim_subset2) - np.nanmin(azim_subset2)) > 90.0:
                            #         print
                            azim_subset3 = azim_subset2.copy()
                            azim_subset2 = unwrap_lons(azim_subset2)
                            inc_subset2 = inc[:,:,look][ilats,ilons]

                            time_resamp = np.nan
                            azim_resamp = np.nan
                            inc_resamp = np.nan

                            for ichan in range(4):
                                tb_subset = tb[:,:,look,ichan][ilats,ilons]
                                ok = np.isfinite(tb_subset)
                                tot_wts = np.sum(wts[ok])
                            
                                if tot_wts > 0.95:
                                    tb_resamp = np.sum((tb_subset*wts)[ok])/tot_wts
                                    lon_resamp = np.sum((lons_subset2*wts)[ok])/tot_wts
                                    if lon_resamp < 0.0:
                                        lon_resamp += 360.0
                                    if lon_resamp > 360.0:
                                        lon_resamp -= 360.0
                                    lat_resamp = np.sum((lats_subset2*wts)[ok])/tot_wts
                                    if ichan == 0:
                                        time_resamp = np.sum((time_subset2*wts)[ok])/tot_wts
                                        azim_resamp = np.sum((azim_subset2*wts)[ok])/tot_wts
                                        if azim_resamp < 0.0:
                                            azim_resamp += 360.0
                                        if azim_resamp > 360.0:
                                            azim_resamp -= 360.0
                                        inc_resamp = np.sum((inc_subset2*wts)[ok])/tot_wts
                                else:
                                    tb_resamp = np.nan
                                    lon_resamp = np.nan
                                    lat_resamp = np.nan
                                    continue

                                ilat_resamp = int(np.round((lat_resamp+90.0)/0.25))
                                ilon_resamp = int(np.round((lon_resamp)/0.25))%1440
                                if ilat_resamp < 0 or ilat_resamp > 720:
                                    print(f'bad lat = {lat_resamp}')
                                    print
                                if (ilat_resamp != ilat_out):
                                    print(f'lat mismatch: ilat_out = {ilat_out}, ilat_resamp = {ilat_resamp}')
                                    print
                                if ilon_resamp < 0 or ilon_resamp > 1440:
                                    print(f'bad lon = {lon_resamp}')
                                    print
                                if ichan == 0:
                                    time_map[ilat_resamp,ilon_resamp,look] = time_resamp
                                    azim_map[ilat_resamp,ilon_resamp,look] = azim_resamp
                                    inc_map[ilat_resamp,ilon_resamp,look] = inc_resamp
                                output_map[ilat_resamp,ilon_resamp,look,ichan] = tb_resamp
            
            #make_plots=True
            if make_plots:
                fig,ax = plot_global_map(stddev_wts[:,:,0],
                                        title='std. dev. wts, fore',
                                        vmin=0,
                                        vmax=0.1,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True) 
                
                fig,ax = plot_global_map(stddev_wts[:,:,1],
                                        title='std. dev. wts, aft',
                                        vmin=0,
                                        vmax=0.1,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True)
                
                fig,ax = plot_global_map(output_map[:,:,0,0],
                                        title='resampled_tbv forward',
                                        vmin=50,
                                        vmax=300,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True) 
                fig,ax = plot_global_map(output_map[:,:,0,1],
                                        title='resampled_tbh forward',
                                        vmin=50,
                                        vmax=300,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True)  
                fig,ax = plot_global_map(output_map[:,:,0,2],
                                        title='resampled_tb3 forward',
                                        vmin=0.0,
                                        vmax=1.0,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True)
                fig,ax = plot_global_map(output_map[:,:,0,3],
                                        title='resampled_tb4 forward',
                                        vmin=0.0,
                                        vmax=1.0,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True)  

                fig,ax = plot_global_map(output_map[:,:,1,0],
                                        title='resampled_tbv aft',
                                        vmin=100,
                                        vmax=300,
                                        cmap='viridis',
                                        plt_colorbar=True,
                                        plt_coastlines=True) 
                plt.show()  
            lat_attrs = coord_attributes_access("latitude", dtype=np.float32)
            lon_attrs = coord_attributes_access("longitude", dtype=np.float32)
            print(f'Assembling {output_file_nc}')
            smap_xr_out = \
                xr.Dataset({'resamp_tbs':(['latitude','longitude','look','pol'],output_map),
                            'resamp_times':(['latitude','longitude','look'],time_map),
                            'resamp_azimuth':(['latitude','longitude','look'],azim_map),
                            'resamp_incidence':(['latitude','longitude','look'],inc_map)},
                            coords={'look':(['look'],np.arange(2)),
                                    'pol':(['pol'],np.arange(4)),
                                    'longitude':(['longitude'],np.arange(1440)*0.25,lon_attrs),
                                    'latitude':(['latitude'],-90.0 + np.arange(721)*0.25,lat_attrs)},
                            attrs=output_attrs,
                        )
            print(f'Writing {output_file_nc}')
            smap_xr_out.to_netcdf(output_file_nc,encoding={ 'resamp_tbs':{'dtype':np.float32,'zlib':True},
                                                            'resamp_times':{'dtype':np.float64,'zlib':True},
                                                            'resamp_azimuth':{'dtype':np.float32,'zlib':True},
                                                            'resamp_incidence':{'dtype':np.float32,'zlib':True}})
            print('.')
            print(f'Wrote {output_file_nc}')
            
        date_to_do += datetime.timedelta(days=1)
    