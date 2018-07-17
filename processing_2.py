"""
The main processing script
----------------------------------------

This script has several processing tasks, enumerated as follows:

1. Read-in all datasets.
2. Find all sampling periods.
3. Allocate all fields for each sampling periods.
4. Write specified fields to output.
5. Call Julia script if cylindrical winds are desired.

To maintain the function-based convention of this computer project. This main script was written as a single function, based from the modules
/:meth:`flightdata` and :meth:`toolbox`. While long functions might seem to deny the actual purpose of function-based scriptting, the use of a main
menu or :meth:`interface` required the use and conection of several functions. For this reason, this and the rest of modules are written based on a single function.

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
from math import radians
# Import glob module to rapidly access all files in a folder.
import glob
# Import scipy module with the packaged of interpolate, in particular, get the function griddata (scipy can be heavy)
from scipy.interpolate import griddata
# Import all functions from toolbox.
from toolbox import *
# Import module to get-in the operating system
import os
# Meteorological python metpy.
import metpy.calc as mpcalc

# Define function.
def read_process_write(filelist,storm):
        """
        ** Name of function says all. **

        *Parameters*

        filelist: `list`
        	Name of storm.
        track: `dict`
        	Dictionary with track. Output of :meth:`flightdata.trackandspeed`.
        storm: `string`
            Name of storm.

        .. note::

            This function makes use of several functions from :meth:`toolbox`, including: :meth:`toolbox.distance`
            It is also important to mention that this function makes use of the programming language Julia by interacting with the operating system and running the Julia script.




        *Returns*
        	temp_axisym.txt:`file` written file with Output fields.

        """
        #Sort filelist.
        filelist=np.sort(filelist)
        print(filelist)
        # Get sampling periods (this will be a dictionary). See the toolbox
        print('Retrieving sampling periods')
        sampleperiods=getsamplingperiods(filelist,2.7)
        sampleperiods={datetime.datetime(1998, 9, 19, 16, 33, 41):datetime.datetime(1998,9,19,21,0,0),datetime.datetime(1998,9,19,21,0,0):datetime.datetime(1998, 9, 20,5, 3,0)}
        omega=7.2921*(10**-5)
        # Iterate over all sampling periods.
        for sampindex,periodskey in enumerate(sampleperiods):

            #Allocate starting (stdt) and ending date (endt). Remeber dt is the convetional short-name for date.
            stdt=periodskey
            endt=sampleperiods[periodskey]

            # Define sampling period string
            if stdt.hour >= 10:
                hh=str(stdt.hour)
            else:
                hh='0'+str(stdt.hour)
            if stdt.day >=10:
                dd=str(stdt.day)
            else:
                dd='0'+str(stdt.day)
            if endt.hour<10:
                fhh='0'+str(endt.hour)
            else:
                fhh=str(endt.hour)
            if endt.day<10:
                fdd='0'+str(endt.day)
            else:
                fdd=str(endt.day)
            period=dd+hh+'-'+fdd+fhh
            print(period)
            # possible user print
            print(stdt,endt)

            # Create new-empty lists.
            lats=[]
            lons=[]

            # Remove outputfile
            os.system('rm temp_axisym.txt')
            dropsincore=0
            print('start filelist loop')
            # Iterate over all files.
            for filename in filelist:



                # Select end-name of file by inspecting filename string. Notice how filename can change how file is read.
                if 'radazm' in filename.split('/')[-1] or 'eol' in filename.split('/')[-1]:
                	end='radazm'
                else:
                	end='avp'
                # Obtain properties of file, i.e., launch time and location into a dictionary (dicc).
                dicc=findproperties(filename,end)

                # Condition to see if current file is in sampling period.
                # Notice how if structure is constructed, condition finds times outside of sampling period and
                # if found outside the sampling period, continue to next file.
                if dicc['Launch Time']<stdt or dicc['Launch Time'] > endt:
                	continue

                # Allocate reading parameters.
                if end =='avp':
                # This section (not used unless indicated by user) is made for the raw datasets or dropsondes downloaded individually.
                    # Read format
                    head=6
                    foot=19
                    # Read of file and allocate into numpy array. with excemption syntax.
                    try:
                        nump=np.genfromtxt(filename,skip_header=head,skip_footer=foot)
                    except:
                        continue
                    # Allocate longitude and latitude.
                    lon=nump[:,11]
                    lat=nump[:,12]
                    # Allocate Temperature (T), Pressure (P), Height (H).
                    T=nump[:,6]
                    P=nump[:,5]
                    H=nump[:,13]
                    #Allocate Relative humidity. (RH)
                    RH=nump[:,7]

                    # Allocate wind magnitude.
                    ur=nump[:,9]
                    # Process wind speed to eliminate false values.
                    ur=cleanu(clean1(ur))

                    # Allocate wind direction.
                    udir=nump[:,8]
                    # Process wind direction.
                    udir=clean1(udir)
                    rfile=nump[:,17]
                    azifile=nump[:,18]
                    # Divide absolute wind speed and wind direction into u-v speeds.
                    u=-ur*np.sin(np.pi*udir/180)
                    v=-ur*np.cos(np.pi*udir/180)

                    # Allocate vertical velocity.
                    w=nump[:,10]


                    # Obtain year/month/day and hour/minute/second values.
                    yymmdd=nump[:,3]
                    hhmmss=nump[:,4]
                    # Convert to single values, hours minutes and seconds.
                    hours,minutes,seconds=timeconversion(hhmmss)

                # Radazm is the label used for the typical dropsonde files used by this computer project.
                # As such, more detailed is given for this part of the conditional argument.
                elif end == 'radazm':
                    # File read parameters.
                    head=16
                    foot=0

                    # Allocate filename fields into a numpy array (nump)
                    nump=np.genfromtxt(filename,skip_header=head,skip_footer=foot)

                    # Allocate Temperature (T), Pressure (P), Height (H).
                    T=nump[:,5]
                    P=nump[:,4]
                    H=nump[:,13]

                    # Allocate time arrays.
                    hours=nump[:,1]
                    minutes=nump[:,2]
                    seconds=nump[:,3]

                    # Allocate Relative Humidity (RH), u and v wind speeds.
                    RH=nump[:,7]
                    u=nump[:,8]
                    v=nump[:,9]

                    # Get wind direction (not used) but skillful user might be interested.
                    udir=nump[:,11]

                    # Obtain vertical windspeed.
                    w=nump[:,-1]
                    rfile=nump[:,17]
                    azifile=nump[:,18]
                    # Get longitude and latitude.
                    lon=nump[:,14]
                    lat=nump[:,15]

                # Clean longitude and latitude from possible nan values.
                lon=clean1(lon)
                lat=clean1(lat)
                # Obtain and round a mean location to proceed to condition to drop environmental dropsondes.
                # and keep inner-core measurements.
                mlon=np.nanmean(lon)
                mlat=np.nanmean(lat)
                lati=np.around(mlat,4)
                longi=np.around(mlon,4)



                # Obtain mean radius (r) and azimuth (theta)
#                r,theta=cart_to_cylindr(mlon,mlat,track,dicc['Launch Time'])

                # Continue condition only if near-inner core dropsonde.
#                if r>290:
#                	continue
                dropsincore+=1
                # print filenames used.
                print(filename)

                # Append longitudes and latitudes to list, useful for quality control on this samplingperiod.
                lats.append(lati)
                lons.append(longi)

                # Clean arrays and save to arrays.
                v_speed=clean2(cleanu(clean1(v)))
                pressure=clean2(cleanu(clean1(P)))
                height=clean2(clean1(H))
                temperature=cleanu(clean2(clean1(T)))
                u_speed=clean2(cleanu(clean1(u)))
                w_speed=cleanu(clean2(clean1(w)))
                relhum=clean2(clean1(RH))
                #if end=='radazm':
	        #        height=getHiso(pressure,temperature,height)
                # Estimate dewpoint using metpy.
                dewpoint=mpcalc.dewpoint_rh((temperature+273)*units.kelvin, relhum/100.)

                # Get equivalent_potential_temperature
                theta_e= equivalent_potential_temperature(pressure*units.hPa, (temperature+273)*units.kelvin, dewpoint)

                # Retrieve potential temperature.
                pot_temp=potential_temperature(pressure*units.mbar,(temperature+273)*units.kelvin)

                # Get storm velocity for this sounding.
                #ustorm,vstorm=stormu(u[1],v[1],dicc['Launch Time'],track[3])

                # Iterate over vectores, we use longitude vector for simplicity but all vectors (height, u_speed, etc) have the same length so for
                # syntax can be used with either vector.
                for j,longi in enumerate(lon):

                    # Get closest date using time arrays. Define datetime object.
                    date=datetime.datetime(dicc['Launch Time'].year,dicc['Launch Time'].month,dicc['Launch Time'].day,int(hours[j]),int(minutes[j]),int(seconds[j]))

                    # Check for nan, if so continue, if observation cannot be located, then all points are useless.
                    if np.isnan(longi) or np.isnan(lat[j]):
                    	continue

                    # Try and obtain radius and azimuth for current observation.
#                    try:
#k                    	r,theta=cart_to_cylindr(longi,lat[j],track,date)
 #                   except:
  #                      continue
                    r=rfile[j]
                    
                    theta=radians(-azifile[j]+90)
            #        print(azifile[j],theta,radians(azifile[j]-90))
                    # Correct u_speed, make it storm relative.
                    u_speed[j]=u_speed[j]#-ustorm
                    v_speed[j]=v_speed[j]#-vstorm

                    # Check for nans in fields. If nans, continue.
                    if np.isnan(r) or np.isnan(theta) or np.isnan(u[j]) or np.isnan(H[j]) or np.isnan(v[j]):
                        continue


                    coriolis_f=(2*omega*np.sin(radians(lat[j])))*10**(2)
                    # After all checks, write file.
                    #print('writing to file')

                    # Open file to append.
                    f=open('temp_axisym.txt','a')

                    # Write all rounded fields to file. While the use of a dictionary and a for loop could synthetize the following line
                    # this approach would use more lines and take more memory (big dictionaries are expensive).
                    f.write(str(np.around(r,5))+'\t'+str(np.around(theta,5))+'\t'+str(np.around(u_speed[j],3))+'\t'+str(np.around(v_speed[j],3))+'\t'+str(np.around(w_speed[j],4))+'\t'
                    +str(np.around(height[j],2))+'\t'+str(np.around(pot_temp[j].magnitude,3))+'\t'+str(np.around(theta_e[j].magnitude,3))+
                    '\t'+str(np.around(pressure[j],3))+'\t'+str(np.around(coriolis_f,4))+'\n')

                    # Close file object (f).
                    f.close()
            print('drops in core '+str(dropsincore))
            print('end of filelist loop')

            # threshold of 6 good dropsondes in sampling period to proceed to call Julia.
            if len(lats)<1:
                os.system('rm temp_axisym.txt')
                continue

            # Why use Julia? you might ask.
            # Why are you not using Julia? I would reply.
            print('Go to Julia @')
            print('julia urutheta.jl')
            os.system('julia urutheta.jl')
            os.system('cp tempjulia.txt outfiles/'+storm+period+'.txt')

            print('starting plotting sequence')
            # Call extra processing or plotting routines.
            #print('python3 3Dfields.py %s %s' % (storm,endt.time()))
#            os.system('python3 gradientwind.py %s %s' % (storm,period))

            # Remove temporary file.
            os.system('rm temp_axisym.txt')
