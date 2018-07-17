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
    ax1.set_ylabel(r'Wind speed [$m$ $s^{-1}$]',fontsize=15)

#    ax1.fill_between([datetime.datetime(2003,9,13,16,0,0),datetime.datetime(2003,9,14,2,0)],[0,100])
    ax1.axvspan(datetime.datetime(2005,9,20,14,14,0),datetime.datetime(2005,9,21,3,0,0),alpha=0.5,color='red')
    ax1.axvspan(datetime.datetime(2005,9,21,15,0,0),datetime.datetime(2005,9,22,1,0,0),alpha=0.5,color='darkslateblue')

   
    ax1.axvspan(datetime.datetime(2005,9,22,14,33,0),datetime.datetime(2005,9,23,3,0),alpha=0.5,color='olive')
    ax1.axvspan(datetime.datetime(2005,9,23,16,0,0),datetime.datetime(2005,9,24,7,0),alpha=0.5,color='green')

    ax1.axvspan(datetime.datetime(2005,9,22,12,0,0),datetime.datetime(2005,9,22,12,25,0),alpha=0.6,color='magenta')
    ax1.axvspan(datetime.datetime(2005,9,23,12,0,0),datetime.datetime(2005,9,23,12,25,0),alpha=0.6,color='magenta')
    ax1.axvspan(datetime.datetime(2005,9,21,10,0,0),datetime.datetime(2005,9,21,10,25,0),alpha=0.6,color='darkorange')
    ax1.axvspan(datetime.datetime(2005,9,22,10,0),datetime.datetime(2005,9,22,10,25,0),alpha=0.6,color='darkorange')
    ax1.plot([datetime.datetime(2005,9,22,12,0,0),datetime.datetime(2005,9,23,12,25,0)],[3,3],'m--',linewidth=3)
    ax1.plot([datetime.datetime(2005,9,21,10,0,0),datetime.datetime(2005,9,22,10,0,0)],[3,3],color='darkorange',linestyle='--',linewidth=3)
    ax1.text(datetime.datetime(2005,9,20,19,30,0),35,"A",fontsize=17,weight='bold',color='white')
    ax1.text(datetime.datetime(2005,9,21,19,0,0),35,"B",fontsize=17,weight='bold',color='white')
    ax1.text(datetime.datetime(2005,9,22,20,0,0),35,"C",fontsize=17,weight='bold',color='white')
    ax1.text(datetime.datetime(2005,9,23,21,0,0),35,"D",fontsize=17,weight='bold',color='white')
    ax1.text(datetime.datetime(2005,9,22,22,0,0),0,"ERC",fontsize=13,weight='bold',color='magenta')
    ax1.text(datetime.datetime(2005,9,21,20,0,0),0,"SEW",fontsize=13,weight='bold',color='darkorange')


    # Put the legend.
    plt.legend(loc='upper left')
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
    plt.legend(loc='upper right')
    plt.xlim([datetime.datetime(2005,9,18,6,0,0),datetime.datetime(2005,9,26,0,0)])
    # Create title
    plt.title(status+' '+storm+' '+year+' intensity history',fontsize=18)
    #Save plot
    plt.savefig('figs/Rita/intensity.png')
    # Show plot
    plt.show()
intensityplot('Rita','2005')
