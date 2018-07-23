"""


This module aims to obtain significant data and meta-data from the flight-level files (downloaded from :cite:`vigh2015flight`) that contain information about the following parameters:

1. Centre track of Tropical Cyclones using the :cite:`willoughby` method with a very good temporal resolution (:math:`\sim 5` minutes).

2. Radius of Maximum Wind derived from wind measurements.

3. Storm speed (wind speed vector of storm motion) :math:`\hat{u_s}=(u_c,v_c)` where :math:`u_c` and :math:`v_c` refer to the x and y components of the storm motion speed vector ( :math:`\hat{u_s}` ).

Simultaneously

This module is briefly divided in the following way:

1. Preamble, loading all packages, functions and modules to be used.
2. Read-in and process of best track-file.
3. Read-in and process of flight-level data.

While several functions were used in the first construction of this module, it now operates based on one single
unified function that carries all the computations and basically does everything. 

.. note::

    This module uses functions from the toolbox:

    :meth:`toolbox.distance`

    And it is the functional basis of most scripts that use this dataset. For example, functions like:

    :meth:`toolbox.getsamplingperiods`, :meth:`toolbox.stormu`, etc.



Main function
==============

"""

# Netcdf module to read-in data
from netCDF4 import Dataset
# Numpy, use of array-objects and mathematical functions.
import numpy as np
# Date and time module to process time easily.
import datetime
# DataFrame processing.
import pandas as pd
# Get distance function from the toolbox.
from toolbox import distance
# Get all global files and matches to a regular expression, also access operating system.
import glob,os
# Trigonometrical functions from math-module.
from math import sin, cos, sqrt, atan2, radians,pi

# Define function
def trackandspeed(storm,year):
    """
    ** Obtain track and speed from flight level data and best track data. **

    *Parameters*

    storm: `string`
    	Name of storm.
    year: `string`
    	Year of storm.

    This function makes use of the function :meth:`toolbox.distance`
    Also, this function was the basis of :ref:`sphx_glr_auto_examples_plot_track.py` so there will be strong similarities.
    Nevertheless, this function comprises more processing and requires further information than that example.

    *Returns*
    	dtime,flat,flon,speeddic:`4-element tuple`

    """
    # Select Flight-level file for our storm and year.
    filen=glob.glob('../Data/'+year+'/'+storm+'/*.nc')
    #Create Dataset netCDF4 Object with file.
    filename=Dataset(filen[0])

    # Read-in best-track file to a Pandas DataFrame
    df=pd.read_csv('track.csv')
    # Select Storm and year.
    df=df[(df['Storm Name']==storm.upper()) & (df['Year']==int(year))]

    # Convert index of dataframe to datetime objects.
    df.index=df['Datetime']
    df.index=pd.to_datetime(df.index)
    #Delete column since it is the index as well.
    del df['Datetime']

    #Create empty lists to be filled in loop.
    btlat=[]
    btlon=[]
    #Create empty dictionary for storm speed motion.
    speeddic={'Datetime':[],'U':[],'V':[]}

    #Date-loop, iteration over all time-steps in best track dataset of our storm.
    #Notice the recurrent use of for loops using an index(i) and the value (dt) of a packed enumerated dataframe.
    for i,dt in enumerate(df.index):
        #Select values for current time-step.
        values=df.loc[dt]
        #Select next time-step to compute speed. (there are 6 h of separation between two values.)
        t1=df.loc[df.index[i+1]]

        # Convert all track longitudes and latitudes to radians.
        lon0=radians(float(values['Longitude'][2:6]))
        lat0=radians(float(values['Latitude'][1:5]))
        lon1=radians(float(t1['Longitude'][2:6]))
        lat1=radians(float(t1['Latitude'][1:5]))

        # Append current latitudes and longitudes to track list.
        btlat.append(float(values['Latitude'][1:5]))
        btlon.append(-float(values['Longitude'][2:6]))

        # Compute difference in longitude and latitude between two points.
        dlon = lon1 - lon0
        dlat = lat1 - lat0

        # Get total distance, x-axis and y-axis distances.
        r=distance(lat0,lon0,lat1,lon1)
        rx=distance(lat0,lon0,lat0,lon1)
        ry=distance(lat0,lon0,lat1,lon0)

        # Get x and y distances to compute speeds.
        x=rx
        y=ry
        # Obtain speed, considering a 6 h separation between two track values.
        u=x/6.
        v=y/6.

        #Current speed is in km/h so conversion to m/s is necessary.
        u=u*1000/3600.
        v=v*1000/3600.

        # Adjust signs of speed to account for x-y +- speeds.
        if dlon > 0:
            u=-u
        if dlat <0:
            v=-v

        # Allocate speed dates and motions in speed-dictionary
        speeddic['Datetime'].append(dt)
        speeddic['U'].append(u)
        speeddic['V'].append(v)

        # Break before reaching end of time-steps in track dataset.
        if i==len(df.index)-2:
            break

    # Possible user print if dictionaries are not your thing.
    #print(speeddic)

    btdate=[]
    # Keys now has all the variable names and output from the flight-level data.
    keys=filename.variables.keys() # print if you are curious.

    # Allocate Radius of Maximum Wind information and dates.
    rmax=np.array(filename['FL_good_radial_leg_flight_level_rmax'])
    rmaxdates=np.array(filename['FL_good_radial_leg_start_Sdatetime'])

    # Datetime index of times where track centre is available.
    dataindex=np.array(filename['FL_WC_wind_center_time_offset'])
    # Latitude and longitude vectors from willoughby-chelmow track. .
    lat=np.array(filename['FL_WC_wind_center_latitude'])
    lon=np.array(filename['FL_WC_wind_center_longitude'])

    # Initial epoch time. (Python standarized date to start "epoch" datetimes)
    t0=datetime.datetime.strptime('1970-01-01 00:00:00','%Y-%m-%d %H:%M:%S')

    # Create flight-level level lists for time (dtime) latitude (flat) and longitude (flon)
    flat=[]
    dtime=[]
    flon=[]

    # Iterate over index of flight-level data.
    for dtindex,cdate in enumerate(dataindex):
        # Obtain date from formatted string (i.e., fromtimestamp)
        date=datetime.datetime.fromtimestamp(cdate).strftime('%Y-%m-%d %H:%M:%S.%f')
        # Second conversion to datetime object.
        date=datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f')

        # Append values to outer lists.
        dtime.append(date)
        flat.append(lat[dtindex])
        flon.append(lon[dtindex])

        # Routine to get speeds, similar algorith, syntax and variable names to
        # those used in the best track case. As such, the following block is un-commented,
        # see above comments for details.
        lon0=radians(flon[dtindex])
        lat0=radians(flat[dtindex])
        lat1=radians(lat[dtindex+1])
        lon1=radians(lon[dtindex+1])
        dlon = lon1 - lon0
        dlat = lat1 - lat0
        r=distance(lat0,lon0,lat1,lon1)
        rx=distance(lat0,lon0,lat0,lon1)
        ry=distance(lat0,lon0,lat1,lon0)
        x=rx*1000
        y=ry*1000

        # Time-step between measurements is variable and not fixed,
        # estimation of time delta is needed. Syntax is below is similar to lines 164-167
        # Allocate next time in iteration.
        newtime=datetime.datetime.fromtimestamp(dataindex[dtindex+1]).strftime('%Y-%m-%d %H:%M:%S.%f')

        # compute delta time or time difference.
        deltatime=datetime.datetime.strptime(newtime, '%Y-%m-%d %H:%M:%S.%f')-date

        # obtain velocity as d/t where time is in seconds.
        u=x/deltatime.seconds
        v=y/deltatime.seconds

        # Correct speed sign for longitude and latitude deltas.
        if dlon > 0:
            u=-u
        if dlat <0:
            v=-v

        # Allocate storm motion results in speed dictionary.
        speeddic['Datetime'].append(dt)
        speeddic['U'].append(u)
        speeddic['V'].append(v)

        # Break condition if appraoching end of timeseries.
        if dtindex==len(dataindex)-2:
            break

    # empty radius of maximum wind dictionary.
    rmwdicc={}

    # loop over rmax datetime objects.
    for i,dt in enumerate(rmaxdates):
        # Error checked allocation of rmws in dictionary
        # Notice how dt is sliced from the first value to the third-to-last values.
        # End of strings is different and can cause errors while reading in a formatted way, as below assumes.
        # Also noticed the error-check syntax is in place to account for possible strings that are poorly formatted in the array.
        try:
            date=datetime.datetime.strptime(dt[0:-3], '%m/%d/%Y %H:%M:%S ')
            # get rmw for current time
            rmwdicc[date]=rmax[i]
        except:
            continue

    #Add rmw results to speed dictonary for compression.
    speeddic['Rmax']=rmwdicc

    #Combine best-track and flight level tracks under one list (i.e. add)
    flat=btlat+flat

    # Return tuple of datetime list, lat and longitudes track list and storm speed dictionary.
    return dtime,flat,flon,speeddic
