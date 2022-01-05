from osgeo import gdal
import numpy as np
import os

#List of numpy arrays we have to process
npdir1=r"\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays"
names1=[r"\Draft_1_1_ft\Draft_1_1_ft_2056",r'\Draft_2_6_ft\Draft_2_6_ft_2056']
Paths=[npdir1+"\\"+name  for name in names1]


#np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223'
#np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\np_arrays'
#np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays'
#np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays\Draft_1_1_ft'
#np_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\Accretion\20211223\np_arrays\Draft_2_6_ft'

ras_dir=r'\\hydro-nas\Team\Projects\5635_SFEI Carbon and GHG\GIS\raster'

#File name
#nam='MTL_2_6_ft_2056'
#Names=['Draft_1_1_ft_2056_years_to_MTL_10']

#Raster template
ras_temp_nam='Layout_Domain.tif'
ras_temp=gdal.Open(os.path.join(ras_dir,ras_temp_nam))
data0=gdal.Open(os.path.join(ras_dir,ras_temp_nam))
geotrans=data0.GetGeoTransform()
proj=data0.GetProjection()
x_pixels,y_pixels=ras_temp.RasterXSize,ras_temp.RasterYSize
driver = gdal.GetDriverByName('GTiff')

for path in Paths:
    nam=path[path.rfind("\\")+1:]
    dst_filename = os.path.join(ras_dir,nam+".tif")
    dataset = None
    dataset = driver.Create(dst_filename,x_pixels, y_pixels, 1,gdal.GDT_Float32)
    #Numpy array
    array=np.load(path+".npy")
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.SetGeoTransform(geotrans)
    dataset.SetProjection(proj)
    dataset.FlushCache()
    dataset=None

