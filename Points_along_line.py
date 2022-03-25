import numpy as np
from shapely.geometry import LineString
from shapely.ops import unary_union
from osgeo import ogr
import geopandas as gpd
import os

#Folfer where we are keeping our vectors
Vector_Path=r"\\hydro-nas\Team\Projects\5625_DSC Climate\GIS\Vector"


distance_delta = 1.64042

##Name of input line shapefile
Input_Name="Line_dummy"
Input = gpd.read_file(os.path.join(Vector_Path,Input_Name+".shp"))
colnames=list(Input.columns)
#colnames.append("Distance")
Points=gpd.GeoDataFrame(columns=colnames)

#Now, let's iterate through all the lines
for line_ind in range(Input.shape[0]):
    line=Input.geometry[line_ind]
    distances = np.arange(0, line.length, distance_delta)
    points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]
    for point_ind in range(len(points)):
        Points_dum = Input.loc[line_ind, Input.columns != "geometry"]
        Points_dum["geometry"] = points[point_ind]
        #Points_dum["Distance"] = distances[point_ind]
        Points = Points.append(Points_dum)
    #multipoint = unary_union(points)
    #Points_dum=Input.loc[line_ind,Input.columns!="geometry"]
    #Points_dum["geometry"]=multipoint


Points.to_file(os.path.join(Vector_Path,Input_Name+"_"+str(round(distance_delta,1))+".shp"))
