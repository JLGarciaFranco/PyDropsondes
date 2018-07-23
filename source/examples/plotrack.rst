.. _sphx_glr_auto_examples_plot_track.py:

Plotting the track of a Tropical Cyclone
=========================================

This example provides a map of the track of a particular tropical cyclone.
The map gives a sense of the longevity of the storm and the strength across its path.

The first and foremost part of this example is the preamble. It is absolutely necessary
that you have all these packages, in particular look out for Basemap!

.. code-block:: python

      # Basemap is a python-module to create accurate and stylish maps.
      from mpl_toolkits.basemap import Basemap
      import matplotlib.pyplot as plt
      # Numpy is our friend, it is the most important module.
      import numpy as np
      # Dates and time module
      import datetime
      # Patches is a secondary tool used to create the handles of the strength-legend.
      import matplotlib.patches as mpatches
      # Pandas, unlike the bear, is used by python to read-in databases, in this case for the track database.
      import pandas as pd
      # We import the function distance from our toolbox.
      from toolbox import distance
      # We import matplotlib, the graphing module tool of python. The abbreviation plt, as np and pd are used as conventional Python docs indicate.
      import matplotlib.pyplot as plt
      # trackhandles is an auxiliary script to create Handle Objects.
      from trackhandles import *


After importing our modules, we make use of the function :ref:`getrack` which reads in the
data and *gets the track* by returning longitude, latitude, windspeed and intensity every 6 h.
This section uses the file *track_w.csv* which must have been provided to you.

.. code-block:: python
      # Get track function of storm and year.
      def getrack(storm,year):
          # we use pandas to read in the data as a Pandas.DataFrame. As pandas-convention dictates, df =dataframe is used as variable name for the object.
          df=pd.read_csv('track_w.csv')
          # We select the storm by also indicating the year. Storm names can be repeated!
          df=df[(df['Storm Name']==storm.upper()) & (df['Year']==int(year))]
          # Pandas DataFrame use indexes, which are useful as they sort our data by date.
          df.index=df['Datetime']
          # Convert index to datetime-objects
          df.index=pd.to_datetime(df.index)
          # empty lists to be filled with values of the track.
          latitudes=[]
          longitudes=[]
          windspeeds=[]
          Intensity_label=[]
          # We iterate over the dates in the pandas dataframe index
          for date in df.index:
              # Select values for current date.
              values=df.loc[date]
              # Append-paste to lists
              latitudes.append(float(values['Latitude'][0:4]))
              longitudes.append(-float(values['Longitude'][1:5]))
              windspeeds.append(float(values['windspeed']))
              Intensity_label.append(values['Intensity'])
          # Return track records.
          return longitudes,latitudes,windspeeds,Intensity_label,df.index




After obtaining the track from the best track datafile, plotting occurs through
the next function :ref:`plotrack`. This function is called by the main script if you decide to plot the track but it can
be used as a standalone if you input the storm and year.



.. code-block:: python
def plotrack(storm,year):
    lon,lat,speedvec,intensity,dates=getrack(storm,year)
    if np.std(lon)>np.std(lat):
        legend_location=(0.5,-0.015)
        columns=3
        plt.figure(figsize=(17,14))
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
            if speed <82:
                dotcolor='yellow'
                size=9
                intensity[index]=r'H$_{1-2}$'

            elif speed >= 82 and speed <112:
                dotcolor='orange'
                size=12
                intensity[index]='H3-4'
            elif speed >= 112:
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
            plt.text(m(lon[index]+.75,lat[index]+1.25)[0],m(lon[index]+.75,lat[index]+1.25)[1],dates[index].day,color='darkgreen',fontsize=10,backgroundcolor='lime')
    m.plot(m(lon,lat)[0],m(lon,lat)[1],'--',color='white')
    plt.title(status+' '+storm+' '+year,fontsize=20)
    c = mpatches.Circle((0.25, 0.25), 0.25, facecolor="orange")
    plt.legend([MObject(),Hurricane3Object(),Hurricane1Object(),TSObject,TDObject,EXObject], ["Major Hurricane (M)","Hurricane Cat. 3-4 (H)","Hurricane Cat. 1-2 (H)","Tropical Storm (TS)","Tropical Depression (TD)","Extratropical Cyclone (EX)"]
        ,handler_map={MObject: MajorHurricaneObjectHandler(),Hurricane3Object: Hurricane3ObjectHandler(),Hurricane1Object: Hurricane1ObjectHandler(),TSObject: TSObjectHandler(),TDObject: TDObjectHandler(),EXObject: EXObjectHandler()},
        title="Legend",bbox_to_anchor=legend_location,ncol=columns)
    plt.show()
