"""
Plot the drift of all dropsondes for a TC
=====================================================

This example gives the typical plot that locates the dropsonde in a storm relative framework.
Typically, this plot is shown in a cartesian coordinate system :math:`(x,y)`, however, this example does it
in cylindrical coordinates :math:`(r,\theta)`, since, actually, to plot the dropsonde location in cartesian coordinates one must first, estimate :math:`r` and :math:`theta`.


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
# Import datetime module for handling date objects.
import datetime
# Import glob module to rapidly access all files in a folder.
import glob
# Import scipy module with the packaged of interpolate, in particular, get the function griddata (scipy can be heavy)
from scipy.interpolate import griddata
# Import all functions from toolbox.
from toolbox import findproperties,getsamplingperiods,cart_to_cylindr,clean1,timeconversion,clean2
# Import module to get-in the operating system
import os

# Define function and arguments (refer to Sphinx)
def plotdrift(filelist,track,storm,show):
        """
        Plot drift function

        """
        # Outpudirectory
        figdir='/home/jlgf/Documents/MRes/Project/figs/'+storm+'/'
        # Create outputdirectory if non existent.
        os.system('mkdir ../figs/'+storm)
        # Get period sampling as a dictionary
        print('Getting sample periods')
        sampleperiods=getsamplingperiods(filelist,2.6)

        print(sampleperiods)
        # First iteration is over the sampling periods to produce one plot for each sampling period.
        # sampleperiods={datetime.datetime(2005, 9, 22, 6, 33, 41):datetime.datetime(2005,9,23,3,0,0),datetime.datetime(2005,9,23,6,0,0):datetime.datetime(2005, 9, 24, 7, 3, 27, 0)}
#        sampleperiods={datetime.datetime(2005, 9, 19, 12, 0, 0):datetime.datetime(2005,9,20,5,0,0),datetime.datetime(2005, 9, 20, 10, 0, 41):datetime.datetime(2005,9,21,3,0,0),datetime.datetime(2005, 9, 22, 5, 33, 41):datetime.datetime(2005,9,23,3,0,0),datetime.datetime(2005,9,23,6,0,0):datetime.datetime(2005, 9, 24, 7, 3, 27, 0)}

        peddict=['A','B','C','D','E','F']
        for sampindex,periodskey in enumerate(sampleperiods):
            # The starting date (sdt) is the key of the dictionary, the value is then the end date (endt)
            pedkey=peddict[sampindex]
            print('Plotting drift in the period :')
            sdt=periodskey
            endt=sampleperiods[periodskey]
            # Print so user knows what time-span corresponds to this plot.
            print(sdt,endt)

            # Create empty lists.
            lats=[]
            lons=[]
            rmaxis=[]
            x=np.array([])
            y=[]
            dates=[]
            dropsondes={}
            maxr=10
            count=0
            #Create figure object with size of 11 megapixels and 9 megapixels.

            # Iterate over dropsondes files.
            for filename in filelist:
# 		possible user print
                print(filename)
                # Establish type of file.
                if 'radazm' in filename.split('/')[-1] or 'eol' in filename.split('/')[-1]:
                	end='radazm'
                else:
                	end='avp'
                # Get information from files.
                dicc=findproperties(filename,end)
                # Evaluate if current file is in sampling period.
                if dicc['Launch Time']>=endt or dicc['Launch Time'] <= sdt:
                # Following print can be enabled to see which files correspond to each plot but this can also saturate printing screen.
                #	print('out of period')
                	continue

                # Read-in data.
                # Notice the following control sequence could have been inserted in line 65-68 above but read-in of data can be computationally expensive so we only read all fields
                # if current file is in sampling period. This control sequence observes the end of th file.
                if end =='avp':

                    # The following description is the same for the end of file =='radazm'
                    # Specify header and footer lengths
                    head=6
                    foot=20

                    # Allocate indexes in file for longitude and latitude.
                    longindex=11
                    latindex=12
                    nump=np.genfromtxt(filename,skip_header=head,skip_footer=foot)

                    # Obtain hours minutes and seconds of measurement.
                    yymmdd=nump[:,3]
                    hhmmss=nump[:,4]
                    hours,minutes,seconds=timeconversion(hhmmss)

                elif end == 'radazm':
                    # as for 'avp'
                    head=16
                    foot=0
                    longindex=14
                    latindex=15
                    nump=np.genfromtxt(filename,skip_header=head,skip_footer=foot)
                    hours=nump[:,1]
                    minutes=nump[:,2]
                    seconds=nump[:,3]

                # Read in file, it is a numpy nd-array, which is suitable for the variable name of nump.

                # Allocate variables.
                lon=nump[:,longindex]
                lat=nump[:,latindex]
                Height=nump[:,13]

                # Clean arrays from possible errors.
                lon=clean2(clean1(lon))
                lat=clean2(clean1(lat))
                Height=clean1(Height)

                # Check for empty arrays, if empty, continue to next file.
                if np.isnan(np.nanmean(lon)) or np.isnan(np.nanmean(lat)):
                	print('NaN error')
                	continue

                # Estimate r and theta from the file.
                r,theta=cart_to_cylindr(np.nanmean(lon),np.nanmean(lat),track,dicc['Launch Time'])

                # If distance is greater than 200 km, this dropsonde is not of interest, then continue to next file.
                if r>150:
                	continue

                # Emtpy lists to allocate for plotting.
                xs=[]
                ys=[]
                rs=[]
                thetas=[]
                # Iteration over longitude array.
                # Break lon in a enumerated tuple where j is the main index.
                for j,longi in enumerate(lon):
                    # Create datetime object using dictionary date and file hours/minutes/seconds
                    date=datetime.datetime(dicc['Launch Time'].year,dicc['Launch Time'].month,dicc['Launch Time'].day,int(hours[j]),int(minutes[j]),int(seconds[j]))

                    # If latitude or longitude are NAN's then continue to next value.
                    if np.isnan(longi) or np.isnan(lat[j]) or np.isnan(Height[j]):
                    	continue

                    # try and get r and theta from values, see sphinx toolbox.
                    try:
                    	r,theta=cart_to_cylindr(longi,lat[j],track,date)
                    except:
                    	continue

                    # Add values to lists.
                    rs.append(r)
                    thetas.append(theta)
                    xs.append(r*np.cos(theta))
                    ys.append(r*np.sin(theta))

                # Condtion to find outer edge of plot and make plot customized to current tiem period.
                # This if only selects the biggest radius (maxr) found across the time=period.
                if len(rs) <1:
                    continue
                if np.nanmax(rs)>maxr:
                    maxr=np.nanmax(rs)
                count+=1
                dropsondes[date]=[thetas,rs]
                # Plot all thetas and radius.
#                ax.plot(thetas,rs,linewidth=3,label=str(date))

            if count <= 10:
                continue
            print('end of filelist loop')

            # Select ticks based on maximum radius.
            rticks=np.arange(0,maxr+10,25)

            plt.figure(figsize=(9,9))
            # Make figure polar
            ax = plt.subplot(111, projection='polar')
            ## Sequence to obtain Radius of Maximum Wind (RMW) from flight-level data.
            rms=track[3]['Rmax']
            ris=0
            counti=0
            # Loop to find all RMW close to this datetime.
            date=dicc['Launch Time']
            for i,key in enumerate(rms):
            	if key>date-datetime.timedelta(hours=3) and key<dicc['Launch Time']+datetime.timedelta(hours=1):
            		ris+=rms[key]
            		counti+=1
            # Average to get a mean RMW of the period.
            if counti!=0:
            	rmax=ris/counti
            else:
                rmax=0
            # Make array of RMW of idntical shape as a plotting array (thetai) to be able to plot RMW.
            rmaxis=[]
            for thetai in np.arange(0,3*np.pi,np.pi/20):
                rmaxis.append(rmax)
            for key in dropsondes.keys():
                thetas,rs=dropsondes[key]
                ax.plot(thetas,rs,linewidth=3)
            # Plot RMW.
            ax.plot(np.arange(0,3*np.pi,np.pi/20.),rmaxis,linewidth=3,color='k')
            print('Dropsondes= '+str(count))
            # Plot settings.
            ax.set_title('Dropsonde drift Hurricane '+storm+' IOP '+str(pedkey),fontsize=16)
            # Set raidus ticks and position.
            ax.set_rticks(rticks)
            ax.set_rlim([0,80])
            ax.set_rlabel_position(135.)

            # Add customized-grid.
            ax.grid(alpha=0.5,linestyle='--')

            #plt.legend()
            # Show Plot.
            plt.savefig('figs/'+storm+'/drift_'+str(periodskey)+'.png')
            if show:
               plt.show()
            plt.close()
