# Notes SMAP resampling

We begin with the publicly available L2C SMAP files.  All processing of SMAP is done in python.  We choose to do this to make the process transparent.  The SMAP data volume is not so large, so the speed concerns with python are not too bad.
 
## Some assumptions that seem to be true:
- The spatial pattern for each SMAP L2C files is the same, but translated in the E-W direction to account for different equator crossing longitudes,
- The footprint geometry, e.g. the azimuth angle of the tilted ellipses is more or less determined by the location in the translated lat/lon grid.

## Approach
For each location in the translated lat/lon grid, we calculate the resampling weights to resample to a 70km circle. Near the equator, we use a 11x11 grid of weights.  As we get nearer the poles, we need more cells in the longitude direction.

Then, for each L2C file, 
- determine the longitude offset
- step through the valid target locations - far enough in from the swath edge so the the resampling works.
- calculate the resampled brightness temperatures and store in the proper location in the output file
- make sure the longitude offset is taken care of.

## storing the weights
For each valid location and look we store:
- The longitude and latitude limits in the L2C source indexing scheme.
- The weights

The weights are stored in separate files for each output latitude. The latitude index, ilat, goes from 20 (-85.0) to 700 (+85.0).



I think maybe 5% of the possible locations are valid, so maybe 50,000 locations, with an average of 800 bytes.  This is 40 MB, so not too big.




