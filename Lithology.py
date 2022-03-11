# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 15:12:55 2019

@author: JDiaz
"""

###First we import the packages. Make sure to have them installed in your Python
import numpy as np
import pandas as pd
import os

#Set parameters manually

if True:
    ###Are we drawing lines this run? True or False
    Lines=True
    
    ###Are we drawing the image?
    Image_draw=False

    #Now we set the size of the size of the page
    #### 17x11 size by default 

    width_size= 1500
    height_size= width_size/1224*792

    ##Now we set the margins in our page

    #Margin in the x left
    mx=100

    #Margin in the x right
    mxr=200

    #Margin in the z
    my=100

    #Margin in the z bottom
    myb=150
    
    ###Significant units in x-axis
    sux=3

    ###Significant units in y-axis
    suy=1



    ###Length of tickmarks

    ltck=10

    ### x-tickmark label offset
    tckxoff=15

    ### y-tickmark label offset
    tckyoff=50



    ###Now we set the width that our wells profiles will have in pixels

    #Well logs width
    wlw=15

    ####The variables xscale and yscale set the size in pixels of the domain where the 
    ####domain of the cross section will be drawn

    #offset of the legend from the right margin
    lgoff=100


    #Size of each pattern box in the legend
    lbz=20
    
    #Distance between legend labels
    delta_leg=10
    
    #Image width
    Iw=200

    #Image height
    Ih=200
    
    #Image vertical position
    Iy=50

    #X-section coordinates
    p0x = 6219784
    p0y = 2190612
    p1x = 6221585
    p1y = 2193854


#Import parameters from CSV

#Current directory
os.chdir(r'P:\Applications\CrossSections')

###Path to spreadsheet
p2s=r"P:\Applications\CrossSections"

if False:
    GraphPar=pd.read_csv(os.path.join(p2s,'GraphicalParameters.csv'),header=None).to_numpy()
    Lines=np.array(GraphPar[0]=='TRUE')[0]
    Image_draw=np.array(GraphPar[1]=='TRUE')[0]
    width_size= int(float(GraphPar[2]))
    height_size= int(float(GraphPar[3]))
    mx=int(float(GraphPar[4]))
    mxr=int(float(GraphPar[5]))
    my=int(float(GraphPar[6]))
    myb=int(float(GraphPar[7]))
    sux=int(float(GraphPar[8]))
    suy=int(float(GraphPar[9]))
    ltck=int(float(GraphPar[10]))
    tckxoff=int(float(GraphPar[11]))
    tckyoff=int(float(GraphPar[12]))
    wlw=int(float(GraphPar[13]))
    lgoff=int(float(GraphPar[14]))
    lbz=int(float(GraphPar[15]))
    delta_leg=int(float(GraphPar[16]))
    Iw=int(float(GraphPar[17]))
    Ih=int(float(GraphPar[18]))
    Iy=int(float(GraphPar[19]))
    p0x =int(float(GraphPar[20]))
    p0y =int(float(GraphPar[21]))
    p1x = int(float(GraphPar[22]))
    p1y =int(float(GraphPar[23]))

#Horizontal scale
xscale=width_size-mx-mxr

#Vertical scale
yscale=height_size-my-myb

#Distance between legend labels
dleg=lbz+delta_leg

#Import wells properties
data = pd.read_csv(os.path.join(p2s,'WellsProps.csv'),encoding='utf-8-sig')
data=data.dropna(axis=1)

##Coordinates of wells
x0=data.filter(items = [0], axis=0)
x0=x0.to_numpy()
x0=x0.reshape(x0.shape[1])
y0=data.filter(items = [1], axis=0)
y0=y0.to_numpy()
y0=y0.reshape(y0.shape[1])


#Layers for each well log

z=pd.read_csv(os.path.join(p2s,'Depths.csv'),encoding='utf-8-sig',header=None)
#Let's drop nas
z=z.dropna(axis=0,how='all')
z=z.dropna(axis=1,how='all')

#Let's transpose
z=z.to_numpy().transpose()



#Geologies
geo=pd.read_csv(os.path.join(p2s,'Lithology.csv'),encoding='utf-8-sig',header=None)
geo=geo.dropna(axis=0,how='all')
geo=geo.dropna(axis=1,how='all')
geo=geo.to_numpy().transpose()


if Lines:
###Import geologic units
    Layers=pd.read_csv(os.path.join(p2s,'GUs.csv'),encoding='utf-8-sig',header=None)
    Layers=Layers.dropna(axis=1,how='all')
    Layers=Layers.iloc[0,:].to_numpy()




    ###Depths of layers
    GUdata=pd.read_csv(os.path.join(p2s,'GUdata.csv'),encoding='utf-8-sig',header=None)
    GUdata=GUdata.dropna(axis=0,how='all')
    GUdata=GUdata.dropna(axis=1,how='all')
    GUdata=GUdata.to_numpy()


###Well names

wellname=list(pd.read_csv(os.path.join(p2s,'WellsProps.csv'),encoding='utf-8-sig').dropna(axis='columns'))



###Well labels horizontal offset
welloff=(np.array([len(i) for i in wellname])*0.3).astype("int")

#Library of lithologies

lith=np.array(['CH','CL','LeveeFill', 'OH', 'OL', 'PT','SC', 'SM', 'SP', 'SW','ML'])

#let´s figure out lithologies in this map

#Flat array with lithologies
lithu=np.array([])
for sl in range(len(geo)):
    lithu=np.append(lithu,geo[sl][geo[sl].astype(str)!=str(np.nan)])

lithu=np.unique(lithu[np.where(lithu!='nan')])

#Patterns library
lithcol=np.array(['url(#CH)','url(#CL)', 'url(#LeveeFill)', 'url(#OH)', 'url(#OL)',\
        'url(#PT)', 'url(#SC)', 'url(#SM)', 'url(#SP)', 'url(#SW)', 'url(#ML)'])

#Transect  coordinates of logs 
x=np.array([])

#Image name

Image='Transect.jpg'

###Image coordinates
Ix=width_size-50-Iw

#Lets get the transect coordinates for the wells
#Lets loop through wells
for well in range(len(x0)):
    qx=x0[well]
    qy=y0[well]
    A=np.matrix([[p1x-p0x,p1y-p0y],
                 [p0y-p1y, p1x-p0x]])
    b=np.matrix([[qx*(p1x-p0x)+qy*(p1y-p0y)],
                 [p0y*(p1x-p0x)-p0x*(p1y-p0y)]])
    A_inv=np.linalg.inv(A)
    xdum=np.matmul(A_inv,b)
    d_dum=int(np.sqrt(xdum[0]**2+xdum[1]**2))
    x=np.append(x,d_dum)

GUdata2=np.array([[],[]])

xdums2=np.empty((int(GUdata.shape[1]/2)-1,GUdata.shape[0]))
xdums2[:]=np.nan
zdums2=np.empty((int(GUdata.shape[1]/2)-1,GUdata.shape[0]))
zdums2[:]=np.nan

if Lines:
####Let's get the transect coordinates for the layers
    for layer in range(len(Layers[Layers!=''])+1):
        xdums=GUdata[~np.isnan(GUdata[:,3*layer]),3*layer]
        ydums=GUdata[~np.isnan(GUdata[:,3*layer+1]),3*layer+1]
        zdums=GUdata[~np.isnan(GUdata[:,3*layer+2]),3*layer+2]
        xdums2_dum=np.array([])
        for point in range(len(xdums)):
            
            qx=xdums[point]
            qy=ydums[point]
            A=np.matrix([[p1x-p0x,p1y-p0y],
                         [p0y-p1y, p1x-p0x]])
            b=np.matrix([[qx*(p1x-p0x)+qy*(p1y-p0y)],
                          [p0y*(p1x-p0x)-p0x*(p1y-p0y)]])
            A_inv=np.linalg.inv(A)
            xdum=np.matmul(A_inv,b)
            d_dum=int(np.sqrt(xdum[0]**2+xdum[1]**2))
            xdums2_dum=np.append(xdums2_dum,d_dum)

        #Let's sort from left to right
        xdums2_dum=xdums2_dum[np.argsort(xdums2_dum)]
        zdums=zdums[np.argsort(xdums2_dum)]
        xdums2[layer][0:len(xdums2_dum)]=xdums2_dum
        zdums2[layer][0:len(zdums)]=zdums
#        GUdata2=np.append(GUdata2,np.array([xdums2,zdums]),axis=1).astype('int')

    #Let's find unique x values of layer points

    x_lay=np.sort(np.unique(xdums2))
    x_lay=x_lay[~np.isnan(x_lay)]

    #Let's create array where we will interpolate and extrapolate everything
    z_lay=np.empty((len(Layers[Layers!=''])+1,len(x_lay)))

    #Let's fill with nas
    z_lay[:]=np.nan

    #Now, we will loop through layers again

    for layer in range(len(Layers[Layers!=''])+1):
        #Let's retrieve the vector of x points that we have for this layer
        x_lay_dum=xdums2[layer]
        x_lay_dum = x_lay_dum[~np.isnan(x_lay_dum)]
        #Let's retrieve the vector of z coordinates that we have for this layer
        z_lay_dum=zdums2[layer]
        z_lay_dum = z_lay_dum[~np.isnan(z_lay_dum)]
        #Now, we'll loop through all the unique x points found in all the layers
        for i in range(len(x_lay)):
            #Now, we have four possible cases
            #The first case is that the x coordinate is to the left of the first point where we have
            #a z coordinate for this layer. In that case, we have to extrapolate to the left, based on the first two
            #points of the layer
            if x_lay[i]<x_lay_dum[0]:
                #Let's get equation of the straight line defined by the first two points of the layer
                fit = np.polyfit(x_lay_dum[:2], z_lay_dum[:2], 1)
                line = np.poly1d(fit)
                z_lay[layer, i]=line(x_lay[i])
            #The second case we have is that the x coordinate is to the right of the last point where we have
            #a z coordinate for this layer. In that case, we have to extrapolate to the right, based on the last two
            #points of the layer
            elif x_lay[i]>x_lay_dum[-1]:
                #Let's get equation of the straight line defined by the last two points of the layer
                fit = np.polyfit(x_lay_dum[-2:], z_lay_dum[-2:], 1)
                line = np.poly1d(fit)
                z_lay[layer, i] = line(x_lay[i])
            #The third case we have is that the x coordinate has a corresponding z coordinate in this layer.
            #In that case, we retrieve the corresponding z value points of the layer
            elif x_lay[i] in x_lay_dum:
                z_lay[layer, i]=np.min(z_lay_dum[x_lay_dum == x_lay[i]])

            #The fourth case we have is that the x coordinate is between two coordinates in this layer.
            #In that case, we interpolate between the values

            else:
                #Index of the point to the left
                ind_left=np.max(np.where(x_lay_dum<x_lay[i]))

                #Index of the point to the right
                ind_right=np.min(np.where(x_lay_dum>x_lay[i]))

                #X coordinates
                x_left=x_lay_dum[ind_left]
                x_right = x_lay_dum[ind_right]

                #Z coordinates
                z_left=z_lay_dum[ind_left]
                z_right = z_lay_dum[ind_right]

                fit = np.polyfit(np.array([x_left,x_right]), np.array([z_left,z_right]), 1)
                line = np.poly1d(fit)
                z_lay[layer, i] = line(x_lay[i])

        #If we are below the top, we have to crop the interpolated values by whichever is lower, the bottom of the
        #upper layer (or surface elevation) or the interpolated value
        if layer>0:
            z_lay[layer]=np.minimum(z_lay[layer],z_lay[layer-1])


    zdums2=z_lay
    xdums2=x_lay




minx=np.sqrt(p0x**2+p0y**2)
maxx=np.sqrt(p1x**2+p1y**2)

#x axis limits

xlim=np.array([minx,maxx])

#Min of z
minz=np.amin(z[np.where(~np.isnan(z))])





maxz=np.amax(z[np.where(~np.isnan(z))])
#z axis limits

if Lines:
    minz=min(minz,np.min(zdums2[~np.isnan(zdums2)]))
    maxz =max(maxz,np.max(zdums2[~np.isnan(zdums2)]))


#ylim=np.array([round(minz/10)*10,round(maxz/10)*10]).astype('int')
ylim=np.array([minz,maxz])

#x tickmarks

###delta x of tickmarks
dtx=10**int(np.log10(xlim[1]-xlim[0])-1)

xticks=np.arange(xlim[0],xlim[1]+dtx,dtx)

xticks0=xticks-minx

#y tickmarks
dty=5
#yticks=(np.arange(ylim[0]-dty,ylim[1]+2*dty,dty)/10).astype("int")*10
yticks=np.arange(int(ylim[0]-dty),int(ylim[1]),dty)
#yticks0=-yticks
yticks0=yticks


#Corrected x coordinates of logs
x2=(x-minx)/(maxx-minx)*xscale+mx-wlw/2



#Corrected x coordinates of tickmarks
xticks=(xticks-minx)/(maxx-minx)*xscale+mx

#width of domain
w=str(int(max(x2)+1.5*wlw+mxr))



#Corrected z coordinates of logs
z2=(z-minz)/(maxz-minz)*yscale+my


if Lines:
    #Corrected z coordinates of layers
    zdums2=(zdums2-minz)/(maxz-minz)*yscale+my
    #Corrected x coordinates
    xdums2=(xdums2-minx)/(maxx-minx)*xscale+mx



#Corrected y coordinates
ylim2=(ylim-minz)/(maxz-minz)*yscale+my
ylim2=ylim2.astype(int)

#Corrected y coordinates of y tickmarks
yticks=(yticks-minz)/(maxz-minz)*yscale+my
yticks=yticks.astype(int)

#x axis limits corrected
xlim2=(xlim-minx)/(maxx-minx)*xscale+mx
xlim2=xlim2.astype(int)
xlim2[1]=xlim2[1]+wlw
xlim2[0]=xlim2[0]-wlw


#height of domain
h=str(int(np.amax(z2[np.where(~np.isnan(z2))])+myb))
z2=int(h)-z2
ylim2=int(h)-ylim2
yticks=int(h)-yticks
if Lines:
    zdums2=int(h)-zdums2

Patterns=open('Patterns.svg',"r")
Patterns=Patterns.read()
f= open('Profile.svg',"w+")
str0='<svg version="1.1"\n     baseProfile="full"\n     width="'+str(width_size)+'" height="'\
+str(height_size)+'"\n     xmlns="http://www.w3.org/2000/svg">\n'

#patterns

str0=str0+Patterns

if Lines:
###Let's draw layers
    for layer in range(len(Layers[Layers!=''])):
        #Top
        #Let's remove nas from the top line of the polygon x and y coordinates
        xdum1=xdums2[~np.isnan(xdums2)]
        zdum1=zdums2[layer][~np.isnan(zdums2[layer])]
        #Let's sort from left to right the y coordinates of the line on top according to their x value for the top line of the polygon
        sort_ind_0 = np.argsort(xdum1)
        zdum1 = zdum1[sort_ind_0]
        # Let's sort x coordinates from left to right for the top line of the polygon
        xdum1 = np.sort(xdum1)

        #Bottom
        # Let's remove nas from the bottom line of the polygon x and y coordinates
        xdum2=xdums2[~np.isnan(xdums2)]
        zdum2 = zdums2[layer + 1][~np.isnan(zdums2[layer+1])]

        #For the bottom line, we sort the points from right to left according to their x coordinate
        sort_ind=np.argsort(xdum2)[::-1]
        zdum2=zdum2[sort_ind]
        xdum2=np.sort(xdum2)[::-1]

        #Now, we format the coordinates for the svg file
        points_dum=np.concatenate((np.array([xdum1,zdum1]),np.array([xdum2,zdum2])),axis=1).astype('int')
        points_dum_str=""
        coldum=lithcol[np.argwhere(lith==Layers[layer])][0][0]
        for i in range(len(points_dum[0])):
#        for i in range(140):
            points_dum_str=points_dum_str+points_dum[0][i].astype('str')+','+points_dum[1][i].astype('str')+' '
    
        str0=str0+'\n  <polygon points="'+points_dum_str+'"'+' fill="'+coldum+'" fill-opacity="0.6" \n'\
                    +'          />\n\n'
    

#Lets loop through wells
for well in range(len(np.array(wellname))):
    
    
    ###well label location
    xwelldum=str(int(x2[well]))
    
    if ~np.isnan(z2[well,0]):
        str0=str0+'\n  <text x="'+str(int(xwelldum)-welloff[well])+'" y="'+str(int(z2[well,0])-5)+\
        '" style="fill:black;stroke:none">'+wellname[well]+'</text>\n'
    
    #Lets scroll through layers
    geodum=geo[well]
    geodum=geodum[np.argwhere(geodum!='')]
    
    for layer in range(len(geodum)):
        
        if ~np.isnan(z2[well,layer+1]):
            #lithology of the layer
            
            litdum=geodum[layer]

            if (len(np.argwhere(lith==litdum))>0)and(~np.isnan(z2[well,layer])):
                #color of the layer
                coldum=lithcol[np.argwhere(lith==litdum)][0][0]
                
                str0=str0+'\n  <rect x="'+xwelldum+'" y="'+str(int(z2[well,layer]))+'" width="'+str(int(wlw))\
                +'" height="'+str(int(z2[well,layer+1]-z2[well,layer]))+'"\n'\
                +'        style="stroke: #000000;\n               fill: '+str(coldum)+';\n'\
                +'              "\n'\
                +'          />\n\n'
    
#Let's plot x axis
yxtick=str(np.amax(yticks))
str0=str0+'\n  <line x1="'+str(xlim2[0])+'" y1="'+str(max(yticks))+'" x2="'+\
str(int(np.amax(xticks)))+'" y2="'+str(max(yticks))+'" stroke="black" />\n'

#Let's plot y axis
xytick=str(xlim2[0])
str0=str0+'\n  <line x1="'+xytick+'" y1="'+str(ylim2[1])+'" x2="'+\
xytick+'" y2="'+str(max(yticks))+'" stroke="black" />\n'

#lets draw x tickmarks
for tick in range(len(xticks)):
    x1dum=str(int(xticks[tick]))
    x2dum=str(int(xticks[tick])-tckxoff)
    y2dum=str(int(yxtick)+ltck)
    y3dum=str(int(yxtick)+2.5*ltck)
    str0=str0+'\n  <line x1="'+x1dum+'" y1="'+yxtick+'" x2="'+\
    x1dum+'" y2="'+y2dum+'" stroke="black" />\n'
    #Label
    str0=str0+'\n  <text x="'+x2dum+'" y="'+y3dum+\
    '" style="fill:black;stroke:none">'+str(xticks0.astype('int')[tick])+'</text>\n' 



#lets draw y tickmarks
for tick in range(len(yticks)):
    x2dum=str(int(xytick)-tckyoff)
    y3dum=str(int(yticks[tick]))
    str0=str0+'\n  <line x1="'+str(int(xytick)-10)+'" y1="'+y3dum+'" x2="'+\
    xytick+'" y2="'+y3dum+'" stroke="black" />\n'
    #Label
    str0=str0+'\n  <text x="'+x2dum+'" y="'+str(int(y3dum)+5)+\
    '" style="fill:black;stroke:none">'+str(yticks0[tick])+'</text>\n'
    
#Let's draw the legend

xleg=str(int(w)-lgoff)
yleg=str(50+Iy+Ih)
for lit in range(len(lithu)):
    #lithology of the layer
    litdum=lithu[lit]
    
    if len(np.argwhere(lith==litdum))>0:    
        #color of the layer
        coldum=lithcol[np.argwhere(lith==litdum)]
        str0=str0+'\n  <rect x="'+xleg+'" y="'+yleg+'" width="'+str(lbz)\
        +'" height="'+str(lbz)+'"\n'\
        +'        style="stroke: #000000;\n               fill: '+str(coldum[0][0])+';\n'\
        +'              "\n'\
        +'          />\n\n'
        
    #Name of the layer
    str0=str0+'\n  <text x="'+str(int(xleg)+lbz+10)+'" y="'+str(int(yleg)+3/4*lbz)+\
    '" style="fill:black;stroke:none">'+lithu[lit]+'</text>\n'    
    yleg=str(int(yleg)+dleg)


###Let's draw image
if Image_draw:    
    str0=str0+'\n  <image href="'+Image+'" x="'+str(Ix)+'" y="'+str(Iy)\
    +'" height="'+str(Ih)+'px" width="'+str(Iw)+'px"/>' 

str0=str0+'</svg>'


line = f.write( str0 )

f.close()
