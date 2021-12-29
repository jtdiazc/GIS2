import numpy as np
from osgeo import gdal
import os

#Rasters directory
ras_dir=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\raster"

#Numpy directory
np_dir=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\np_arrays"

layer_name='SLR_2_6'


ds = gdal.Open(os.path.join(ras_dir,layer_name+'.tif'))
myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
#myarray=np.flip(myarray, 0)
ds=None
#na_value=np.sort(np.unique(myarray))[0]
#myarray[myarray==na_value]=np.nan

np.save(os.path.join(np_dir,layer_name+'.npy'),myarray)