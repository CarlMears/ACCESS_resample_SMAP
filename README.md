# SMAP resampling

## This project contains python code to resample L2C SMAP radiance data to the ACCESS grid and footprints

The starting point is the RSS-produced SMAP L2C orbital files.  In these files, the TOA radiances are already resampled onto an 0.25-degree X 0.25-degree Earth-fixed grid.  The footprints are the native SMAP size (47 x 39 km), and are canted at the native Earth Azimuth Angle. The L2C grid is a "cell centered" grid, with grid centers offsets by 0.125 degrees from the 0.25 degree ACCESS grid.

The ACCESS footprints are centered at regular 0.25 degree locations, i.e. 0.00, 0.25, 0.50, etc.  The ACCESS footprints for SMAP are 70 km gaussians to correspond to the ACCESS footprints for AMSR2 and SSM/I.

## Approach

- The spatial pattern for each of the SMAP L2C files is the same, but translated in the E-W direction to account for different equator crossing longitudes.
- The footprint geometry, e.g. the azimuth angle of the tilted ellipses is determined by the location in the translated lat/lon grid.
- We calculate the resampling weight on the L2C grid, taking into account the azimuth angle of each source footprint.  Note that this means that the weights are different for the fore and aft looks because of the azimuth angle differences.

## Calculating the Weights
- In the SMAP L2C files, there is some amount of jitter in the azimuth angle for individual resampled footprints, and there are regions of missing data in each orbit.  To get a representative map of azimuth angles, we average over several orbits in the jupyter notebook: 
    ~~~
    calc_average_location_error_stats_over_multiple_orbits.ipynb
    ~~~

- We then calculate the backus-gilbert weights for each output location.  The number of source footprints that are considered in the resampling calculations depends on secant(latitude).  At the equator, a 10 x 10 array of source footprints is considers.  At 85.0 degrees (North or South) we consider a 60 x 10 array of source footprints.  At the equator, the smoothing factor is set to a small number (0.000001) so that only a minimal amount of smoothing is done.  Poleward of 75 degrees, this needs to be increased by a factor of 10 to ensure numerical stability because of the increasing overlap of the source footprints.  The calculations in done in the code:
    ~~~
    calc_resample_wts_SMAP_to_grid_circular_footprint.py
    ~~~

- There are separate output files for each ilat.  The files are named output files are named with names like:
    ~~~
    L2C_weights_692.70km.v2.nc,
    ~~~
    where the "692" is the ilat value. These files a currently available to the public by request to support@remss.com.

## Processing SMAP orbits

Processing of the SMAP orbits is performed by
   ~~~ 
   resample_SMAP_L2C_files.py
   ~~~
Several named arguments are required by the routine:
* --smap_wt_root: location of the weight files
* --output_root: location of the output files
* --start_date: (in yyyy-mm-dd) 
* --end_date: (in yyy-mm-dd) 
* --verbose: Include to generate more descriptive output to the terminal

There are several ".bat" files in the repo with examples of calling the routine.