from osgeo import gdal
import numpy as np
import os

#Directory with numpy arrays

np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays\Draft_1_1_ft'

ras_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\raster'

#File name
nam='Draft_1_1_ft_2018'

#Raster template
ras_temp_nam='Layout_Domain.tif'
ras_temp=gdal.Open(os.path.join(ras_dir,ras_temp_nam))
dst_filename = os.path.join(ras_dir,nam+".tif")
x_pixels,y_pixels=ras_temp.RasterXSize,ras_temp.RasterYSize
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(dst_filename,x_pixels, y_pixels, 1,gdal.GDT_Float32)

#Numpy array
array=np.load(os.path.join(np_dir,nam+".npy"))

dataset.GetRasterBand(1).WriteArray(array)

data0=gdal.Open(os.path.join(ras_dir,'Layout_Domain.tif'))

geotrans=data0.GetGeoTransform()

proj=data0.GetProjection()

dataset.SetGeoTransform(geotrans)

dataset.SetProjection(proj)

dataset.FlushCache()

dataset=None

