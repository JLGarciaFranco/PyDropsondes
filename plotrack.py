# -*- coding: utf-8 -*-
"""
Plotting example 1
========================

The gallery is capable of transforming Python files into reStructuredText files
with a notebook structure. For this to be used you need to respect some syntax
rules.

It makes a lot of sense to contrast this output rst file with the
:download:`original Python script <plot_notebook.py>` to get better feeling of
the necessary file structure.

Anything before the Python script docstring is ignored by sphinx-gallery and
will not appear in the rst file, nor will it be executed.
This Python docstring requires an reStructuredText title to name the file and
correctly build the reference links.

Once you close the docstring you would be writing Python code. This code gets
executed by sphinx gallery shows the plots and attaches the generating code.
Nevertheless you can break your code into blocks and give the rendered file
a notebook style. In this case you have to include a code comment breaker
a line of at least 20 hashes and then every comment start with the a new hash.

As in this example we start by first writing this module
style docstring, then for the first code block we write the example file author
and script license continued by the import modules instructions.
"""

try:
	from mpl_toolkits.basemap import Basemap
except ImportError:
	print("It would seem that you do not have Basemap properly installed. Example failed. ")
	quit()
import matplotlib.pyplot as plt
# setup Lambert Conformal basemap.
# set resolution=None to skip processing of boundary datasets.
from netCDF4 import Dataset
import numpy as np
import datetime
import matplotlib.patches as mpatches
import pandas as pd
from toolbox import distance
import glob,os
from math import sin, cos, sqrt, atan2, radians,pi
import matplotlib.pyplot as plt
from trackhandles import *
def getrack(storm,year):
    df=pd.read_csv('track_w.csv')
    df=df[(df['Storm Name']==storm.upper()) & (df['Year']==int(year))]
    df.index=df['Datetime']
    df.index=pd.to_datetime(df.index)
    del df['Datetime']
    btlat=[]
    btlon=[]
    windspeed=[]
    Intensity_label=[]
    speeddic={'Datetime':[],'U':[],'V':[]}
    for i,dt in enumerate(df.index):
        values=df.loc[dt]
        btlat.append(float(values['Latitude'][0:4]))
        btlon.append(-float(values['Longitude'][1:5]))
        windspeed.append(float(values['windspeed']))
        Intensity_label.append(values['Intensity'])
    return btlon,btlat,windspeed,Intensity_label,df.index
def plotrack(storm,year):
    lon,lat,speedvec,intensity,dates=getrack(storm,year)
    if np.std(lon)>np.std(lat):
        legend_location=(0.869,0.759)
        columns=1
        plt.figure(figsize=(14,10))
    else:
        legend_location=(1.05, 1)
        columns=1
        plt.figure(figsize=(14,12))
    m = Basemap(llcrnrlon=np.nanmin(lon)-5,llcrnrlat=np.nanmin(lat)-2.5,urcrnrlon=np.nanmax(lon)+2.*np.std(lon),urcrnrlat=np.nanmax(lat)+2.5,
                projection='lcc',lon_0=np.mean(lon),lat_0=np.mean(lat),
                resolution ='l',area_thresh=1000.)
    m.bluemarble()
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    parallels = np.arange(0.,81,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False],color='wheat')
    meridians = np.arange(0.,351.,15.)
    m.drawmeridians(meridians,labels=[True,False,False,True],color='white')
    if len(speedvec)>22:
        divmodulator=4
    if len(speedvec)>14:
        divmodulator=3
    else:
        divmodulator=2
    status=''
    daylist=[]
    for index,speed in enumerate(speedvec):
        if speed > 34 and speed < 64:
            dotcolor='aqua'
            size=6.5
            if status!="Hurricane":
                status="Tropical Storm"
        elif speed >=64:
            status='Hurricane'
            if speed <96:
                dotcolor='yellow'
                size=9
                intensity[index]=r'H$_{1-2}$'

            elif speed >= 96 and speed <136:
                dotcolor='orange'
                size=12
                intensity[index]='H3-4'
            elif speed >= 136:
                dotcolor='red'
                size=18
                intensity[index]='M'
        else:
            if status!="Hurricane" and status!="Tropical Storm":
                status="Depression"
            dotcolor='blue'
            size=3.5
        if intensity[index]==' EX':
            dotcolor='magenta'
            size=5
        m.scatter(m(lon[index],lat[index])[0],m(lon[index],lat[index])[1],s=size*7.5,c=dotcolor)
            #plt.text(m(lon[index],lat[index])[0],m(lon[index],lat[index])[1],intensity[index],color='white',fontsize=12)
        if dates[index].day not in daylist:
            daylist.append(dates[index].day)
            plt.text(m(lon[index]+.75,lat[index]+1.)[0],m(lon[index]+.75,lat[index]+1.)[1],dates[index].day,color='white',fontsize=14,backgroundcolor='black')
    m.plot(m(lon,lat)[0],m(lon,lat)[1],'--',color='white')
    plt.title(status+' '+storm+' '+year,fontsize=20)
    c = mpatches.Circle((0.25, 0.25), 0.25, facecolor="orange")
    plt.legend([MObject(),Hurricane3Object(),Hurricane1Object(),TSObject,TDObject,EXObject], ["Major Hurricane (M)","Hurricane Cat. 3-4 (H)","Hurricane Cat. 1-2 (H)","Tropical Storm (TS)","Tropical Depression (TD)","Extratropical Cyclone (EX)"]
        ,handler_map={MObject: MajorHurricaneObjectHandler(),Hurricane3Object: Hurricane3ObjectHandler(),Hurricane1Object: Hurricane1ObjectHandler(),TSObject: TSObjectHandler(),TDObject: TDObjectHandler(),EXObject: EXObjectHandler()},
        title="Legend",loc='upper right')#,bbox_to_anchor=legend_location,ncol=columns)
    plt.show()
#plotrack('Isabel','2003')
