import numpy as np
import xarray as xr
from pathlib import Path 
import os
import glob
import datetime

if os.name == 'nt':
    root = Path("//athena/public/ftp/smap/SSS/V05.0/FINAL/L2C") 
elif os.name == 'posix':
    root = Path("/mnt/athena/public/ftp/smap/SSS/V05.0/FINAL/L2C")

def find_SMAP_L1C_filenames_for_orbit(*,
                                      date: int,
                                      path : Path = root):
    # date is a datetime.date object
    # path is a pathlib.Path object
    # returns a Path object
    year = date.year
    month = date.month
    day = date.day
    path_to_file = path / f'{year}/{month:02d}'
    date_string = f'{year:04d}{month:02d}{day:02d}'

    files = glob.glob(f'{path_to_file}/RSS_SMAP_SSS_L2C_r*_{date_string}T*_FNL_V05.0.nc')
    paths = [Path(file) for file in files]
    return paths

def read_SMAP_L1C_file(*,
                       path: Path):
    # path is a pathlib.Path object
    # returns an xarray.Dataset object
    ds = xr.open_dataset(path,decode_times=False)
    return ds





