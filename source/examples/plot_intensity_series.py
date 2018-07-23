"""
Plotting the intensity metrics of a Tropical Cyclone
=====================================================

This example gives the typical plot of pressure and wind speed from the best track dataset.

Specifically, the `best track <https://www.nhc.noaa.gov/data/#hurdat>`_ dataset, provides the maximum sustained winds at 10 m altitude (hereafter :math:`U_{10}`)
and the minimum surface pressure (hereafter :math:`P_{min}`).

"""
# matplotlib is the main module in python for plotting.
import matplotlib.pyplot as plt
# Numpy is our friend, it is the most important module.
import numpy as np
# Dates and time module
import datetime
# Pandas, unlike the bear, is used by python to read-in databases, in this case for the track database.
import pandas as pd
# Get track function of storm and year.
def intensityplot(storm,year):
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
    pressures=[]
    Intensity_label=[]
    # Iterate over the dates in the pandas dataframe index
    for date in df.index:
        # Select values for current date.
        values=df.loc[date]
        # Append-paste to lists
        latitudes.append(float(values['Latitude'][0:4]))
        longitudes.append(-float(values['Longitude'][1:5]))
        windspeeds.append(float(values['windspeed']))
        pressures.append(float(values['pressure']))
        Intensity_label.append(values['Intensity'])
    if np.nanmax(windspeeds)>64:
        status='Hurricane'
    elif np.nanmax(windspeeds)>34:
        status='Tropical Storm'
    else:
        status='Tropical Depression'
    windspeeds=np.array(windspeeds)
    windspeeds=windspeeds*0.514444444

    # Create the figure, with a certain size
    fig, ax1 = plt.subplots(figsize=(12,7))
    # Plot wind speed.
    ax1.plot(df.index,windspeeds,'k',label=r'$U_{10}$',linewidth=3)
    # Wind label.
    ax1.set_ylabel(r'Wind speed [$m s^{-1}$]',fontsize=15)
    # Put the legend.
    plt.legend()
    # Make the twin-axis
    ax2=ax1.twinx()
    # Plot the pressure field
    ax2.plot(df.index,pressures,'b--',label=r'$P_{min}$',linewidth=3)
    # Put pressure label
    ax2.set_ylabel('Pressure (mb)',fontsize=15)
    # X-axis label is time
    ax1.set_xlabel('Date',fontsize=15)
    #Put grid-on and the legend.
    plt.grid()
    plt.legend()
    # Create title
    plt.title(status+' '+storm+' '+year,fontsize=18)
    #Save plot
    plt.savefig('figs/Ivan04intensity.png')
    # Show plot
    plt.show()
intensityplot('Ivan','2004')
