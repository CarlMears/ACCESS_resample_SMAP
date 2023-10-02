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
from collections import OrderedDict

def flatten_4d_to_2d(z):
    z = z.reshape(z.shape[0],z.shape[1],-1)
    z = z.reshape(-1,z.shape[-1])
    return z

def flatten_2d_to_1d(z):
    z = z.reshape(-1)
    return z



def find_nearby_footprint(*,target_lat,
                            target_lon,
                            source_lats,
                            source_lons,
                            delta_ilat,
                            verbose = False):

    # average the source lons along the NS diretion

    mn_source_lons = np.nanmean(source_lons,axis=0)

    # find the closest source lon to the target lon
    dlon = (180.0 + mn_source_lons - target_lon) % 360.0 - 180.0
    abs_dlon = np.abs(dlon)
    abs_dlon[np.isfinite(dlon) == False] = 9999.0
    ok = np.where(abs_dlon < 0.2)
    if len(ok[0]) == 0:
        return (None,None,None,None)
    min_dlon_index = ok[0][0]
    if verbose:
        print(f'min_dlon_index: {min_dlon_index}')
    min_dlon = abs_dlon[min_dlon_index]

    if min_dlon_index < 60.0:
        min_dlon_index += 1440
        min_dlon2 = abs_dlon[min_dlon_index]
        if verbose:
            print(f'min_dlon2: {min_dlon2}')
    if min_dlon_index > 1500:
        min_dlon_index -= 1440
        min_dlon2 = abs_dlon[min_dlon_index]
        if verbose:
            print(f'min_dlon2: {min_dlon2}')


    if verbose:
        print(f'min_dlon: {min_dlon}')
        print(f'min_dlon_index: {min_dlon_index}')
    
    # find the closest source lat to the target lat
    lats_at_min_dlon = source_lats[:,min_dlon_index]
    dlat = np.abs(lats_at_min_dlon - target_lat)
    dlat[np.isfinite(dlat) == False] = 9999.0
    ok = np.where(dlat < 0.2)
    if len(ok[0]) == 0:
        return (None,None,None,None)
    
    min_dlat_index = ok[0][0]
    min_dlat = dlat[min_dlat_index]
    if verbose:
        print(f'min_dlat: {min_dlat}')
        print(f'min_dlat_index: {min_dlat_index}')

    delta_ilon = int(np.round(delta_ilat/np.cos(np.deg2rad(target_lat))))
                     
    lon_index_range = (min_dlon_index-delta_ilon,min_dlon_index+delta_ilon+1)
    lat_index_range = (min_dlat_index-delta_ilat,min_dlat_index+delta_ilat+1)

    assert lon_index_range[0] >= 0
    assert lat_index_range[0] >= 0
    assert lon_index_range[1] < source_lons.shape[1]
    assert lat_index_range[1] < source_lons.shape[0]

    lon_indices = np.arange(lon_index_range[0],lon_index_range[1])
    lat_indices = np.arange(lat_index_range[0],lat_index_range[1])

    lon_indices,lat_indices = np.meshgrid(lon_indices,lat_indices)

    return (lon_indices,lat_indices,min_dlon,min_dlat)

if __name__ == "__main__":

    verbose = True
    smoothing_factor = 0.00001
    delta_ilat = 5
    output_path = Path('M:/job_access/resampling/SMAP/weights_for_L2C_files_v3')

    #load the previously averaged values of the azimuth angle
    da = xr.open_dataarray(output_path / 'mean_eaa_map.nc')
    mean_eaa_map = da.data
    input_lats= np.arange(-90.0+0.125,89.9,0.25)
    input_lons = np.arange(0.125-15.0,359.99+15.0,0.25)
    input_lat_grid,input_lon_grid = np.meshgrid(input_lats,input_lons,indexing='ij')

    for itar_lat in range(55,54,-1):
        g_ordered_dict_0 = OrderedDict()
        g_ordered_dict_1 = OrderedDict()
        g_ordered_dict = [g_ordered_dict_0,g_ordered_dict_1]

        target_lat = -90.0 + itar_lat/4.0
        delta_ilon = int(np.round(delta_ilat/np.cos(np.deg2rad(target_lat))))
        if delta_ilon > 30:
            delta_ilon=30
        num_x = 2*delta_ilon 
        num_y = 2*delta_ilat 
        num_output_lat = 1440
        num_looks = 2

        #zero the output arrays.  
        itar_lat_array = np.full(num_output_lat,-999,dtype=np.int32)
        itar_lon_array = np.full(num_output_lat,-999,dtype=np.int32)
        wts_array = np.full((num_output_lat,num_looks,num_y,num_x),np.nan,dtype=np.float64)
        lat_indices_array = np.full((num_output_lat,num_y,num_x),-999,dtype=np.int32)
        lon_indices_array = np.full((num_output_lat,num_y,num_x),-999,dtype=np.int32)
        num_valid_locations = 0

        for itar_lon in range(0,1440):

            target_lon = itar_lon/4.0
            #print(f'target_lat: {target_lat}, target_lon: {target_lon}')

            latitude_center = input_lat_grid[itar_lat,itar_lon]
            longitude_center = input_lon_grid[itar_lat,itar_lon+60]
            mean_eaa = mean_eaa_map[itar_lat,itar_lon+60,0]

            if np.isfinite(mean_eaa) == False:
                continue
            #print(f'itar_lat: {itar_lat}, itar_lon: {itar_lon}, latitude_center: {latitude_center:.3f}, longitude_center: {longitude_center:.3f}')


            range_x = np.arange(-delta_ilon+itar_lon+60,delta_ilon+itar_lon+60)
            #print(range_x)

            range_y = np.arange(-delta_ilat+itar_lat,delta_ilat+itar_lat)
            #print(range_y)

            x_indices,y_indices = np.meshgrid(range_x,range_y)
            lat_indices_array[num_valid_locations,:,:] = y_indices
            lon_indices_array[num_valid_locations,:,:] = x_indices

            stored_valid_wts = False
            for look in range(2):
                near_eaa = mean_eaa_map[y_indices,x_indices,look]
                near_lons = input_lon_grid[y_indices,x_indices]
                near_lats = input_lat_grid[y_indices,x_indices]
                

                if np.sum(np.isfinite(near_eaa)) < near_eaa.size:
                    #if verbose:
                        #print('NaNs in near_eaa')  
                    continue
        
                grid = LocalEarthGrid(center_lon = target_lon,center_lat = target_lat,
                                    delta_x = 1.,delta_y = 1.,
                                    num_x = 1801,num_y = 1201,
                                    name='Single',
                                    do_transform=False)
        
                x0,y0 = grid.projection(near_lons,near_lats)
                x0 = x0/1000.0
                y0 = y0/1000.0

                numy = x0.shape[0]
                numx = x0.shape[1]

                # These are the array with the info about the self overlaps
                g_array = np.full((numy,numx,numy,numx),0.0)
                ratio_array = np.full((numy,numx,numy,numx),np.nan)
                distance_array = np.full((numy,numx,numy,numx),np.nan)
                num_nonzero = 0
                num_found = 0
                num_computed = 0

                # These are the arrays with the info about the overlap with the target footprint
                v_array = np.full((numy,numx),0.0)

                g_target = Gaussian(0.0,0.0,70.0,70.0,0.0)
                verbose = False
                for ix1 in range(numx):
                    for iy1 in range(numy):
                        g1 = Gaussian(x0[iy1,ix1],y0[iy1,ix1],47.0,39.0,float(near_eaa[iy1,ix1]))
                        self_overlap = g1.overlap(g1)
                        v_array[iy1,ix1] = g1.overlap(g_target)
                        if verbose:
                            print(f'ix1: {ix1}, iy1: {iy1}, {v_array[iy1,ix1]:.8f}')

                        for ix2 in range(numx):
                            for iy2 in range(numy):
                                distance = np.sqrt(np.square(x0[iy1,ix1] - x0[iy2,ix2]) + np.square(y0[iy1,ix1] - y0[iy2,ix2]))
                                if distance < 120.0:
                                    key = f'{lat_indices_array[num_valid_locations,iy1,ix1]:04d}'
                                    key = key + f'-{lon_indices_array[num_valid_locations,iy1,ix1]:04d}'
                                    key = key + f'-{lat_indices_array[num_valid_locations,iy2,ix2]:04d}'
                                    key = key + f'-{lon_indices_array[num_valid_locations,iy2,ix2]:04d}'
                                    key = key + f'-{look:02d}'
                                    try:
                                        overlap = g_ordered_dict[look][key]
                                        num_found += 1
                                    except KeyError:
                                        g2 = Gaussian(x0[iy2,ix2],y0[iy2,ix2],47,39,float(near_eaa[iy2,ix2]))
                                        overlap = g1.overlap(g2)
                                        g_ordered_dict[look][key] = overlap
                                        num_computed += 1
                                    ratio = overlap/self_overlap

                                    if ratio > 0.0001:
                                        g_array[iy1,ix1,iy2,ix2] = overlap
                                        ratio_array[iy1,ix1,iy2,ix2] = ratio
                                        distance_array[iy1,ix1,iy2,ix2] = distance
                                        num_nonzero += 1
                                    if verbose:
                                        print(f'ix1: {ix1}, iy1: {iy1}, ix2: {ix2}, iy2: {iy2}, {overlap:.8f}, distance = {distance:.2f}, ratio = {overlap/self_overlap:.5f}')
                print(f'Num Computed: {num_computed}, Num Found: {num_found}, Size of Dict: {len(g_ordered_dict[look])}')
                print
                
                #flatten out the x and y dimensions so the calculations can be done with matrices
                g_array = flatten_4d_to_2d(g_array).astype(np.float64)
                ratio_array = flatten_4d_to_2d(ratio_array)
                distance_array = flatten_4d_to_2d(distance_array)
                v_array = flatten_2d_to_1d(v_array).astype(np.float64)

                u = np.ones(g_array.shape[0],dtype=np.float64)

                # add a smoothing factor to the diagonal elements
                # this is the capital V in the equations
                eye = np.eye(g_array.shape[0],dtype=np.float64)
                g_array = g_array + smoothing_factor*eye

                g_inv = np.linalg.inv(g_array)

                g_inv_v = np.matmul(g_inv,v_array)
                u_g_inv_v = np.matmul(u,g_inv_v)
                numerator = 1.0 - u_g_inv_v

                g_inv_u = np.matmul(g_inv,u)
                demoninator = np.matmul(u,g_inv_u)

                wts = np.matmul(g_inv,(v_array + (numerator/demoninator)*u))

                wts = wts.reshape(numy,numx)

                if num_valid_locations == 0:
                    plt.imshow(wts)

                # write the weights to a file
                itar_lat_array[num_valid_locations] = itar_lat
                itar_lon_array[num_valid_locations] = itar_lon
                wts_array[num_valid_locations,look,:,:] = wts
                print(f'itar_lat: {itar_lat}, itar_lon: {itar_lon}, look: {look}, num_valid_locations: {num_valid_locations},total_wts: {np.sum(wts)}')
                print(f'std.dev. wts = {np.std(wts)}')
                stored_valid_wts = True
                print
            if stored_valid_wts:
                num_valid_locations += 1 

        if num_valid_locations > 0:
            itar_lat_array = itar_lat_array[0:num_valid_locations]
            itar_lon_array = itar_lon_array[0:num_valid_locations]
            wts_array = wts_array[0:num_valid_locations,:,:,:]
            lat_indices_array = lat_indices_array[0:num_valid_locations,:,:]
            lon_indices_array = lon_indices_array[0:num_valid_locations,:,:]

            print(f'itar_lat_array.shape: {itar_lat_array.shape}')
            print(f'itar_lon_array.shape: {itar_lon_array.shape}')
            print(f'wts_array.shape: {wts_array.shape}')
            out_xr = xr.Dataset({'weights': (['locs','look','y','x'],wts_array),
                                 'lat_indices': (['locs','y','x'],lat_indices_array),
                                 'lon_indices': (['locs','y','x'],lon_indices_array)},
                                coords={'lat': (['lat'], itar_lat_array),
                                        'lon': (['lon'], itar_lon_array),
                                        'look': (['look'], np.arange(2)),
                                        'y': (['y'], np.arange(-delta_ilat+1,delta_ilat+1)),
                                        'x': (['x'], np.arange(-delta_ilon+1,delta_ilon+1)),
                                        'locs': (['locs'], np.arange(num_valid_locations))},
                                attrs={'description': 'weights for resampling SMAP L1C to 70km grid',
                                       'smoothing_factor': f'{smoothing_factor}'})
            out_xr.to_netcdf(output_path / f'L2C_weights_{itar_lat:03d}.70km.v2.nc')

        print




        



   