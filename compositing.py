"""
Compositing TCs
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
# Import glob module to rapidly access all files in a folder.
import glob
# Import scipy module with the packaged of interpolate, in particular, get the function griddata (scipy can be heavy)
from scipy.interpolate import griddata

from flightdata import trackandspeed

# Import all functions from toolbox.
from toolbox import *
# Import module to get-in the operating system
import os
# Meteorological python metpy.
import metpy.calc as mpcalc

def get_Rmax(track,startdate,endate):
    ## Sequence to obtain Radius of Maximum Wind (RMW) from flight-level data.
    rms=track[3]['Rmax']
    ris=[]
    counti=0
    rms=dict(rms)
    # Loop to find all RMW close to this datetime.
    #for key in period.keys():
    #    startdate=key
    #    endate=period[key]
    for i,key in enumerate(rms):
    	if key>startdate and key<endate:
    		ris.append(rms[key])
    return np.nanmean(ris)

directorynames=['/gps.qc.eol/GIV/','/ublox.qc.eol/GIV/','/gps.qc.eol/P-3.43/','/ublox.qc.eol/P-3.43/','/gps.qc.eol/P-3.42/','/ublox.qc.eol/noaa.P-3/','/ublox.qc.eol/P-3.42/','/gps.qc.eol/USAF/','/ublox.qc.eol/USAF/']
home=os.getcwd()
home=home+'/'
datadir=home+'/'
for i in range(1999,2013):

    year=str(i)
    yrdir=home+'/Data/'+year+'/'
    stdirs=os.listdir(yrdir)

    for storm in stdirs:

        filelist=[]

        for direct in directorynames:

            filelist=filelist+glob.glob(home+'/Data/'+year+'/'+storm+direct+'*')
            filelist=glob.glob(home+'/Data/'+year+'/Hawk/*')+filelist
        filelist=np.sort(filelist)
        print(storm,year)
        if storm=='Hawk':
            continue
        sampleperiods=getsamplingperiods(filelist,3.)
        try:
            os.system('cp /media/jlgf/Seagate\ Expansion\ Drive/FlightData/v1.1/*'+year+'*'+storm.upper()+'* '+'Data/'+year+'/'+storm+'/')
            track=trackandspeed(storm,year)
        except:
            print('No FlightData for '+storm+' '+year)
            continue


        os.system('rm '+'Data/'+year+'/'+storm+'/*.nc')
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
            inten,Ic=periodI(storm,stdt.year,stdt,endt)
            if np.isnan(inten) or np.isnan(np.mean(Ic)):
                continue
            line=storm+','+str(stdt)+','+str(endt)+','+str(inten)+','+str(Ic[0])+'\n'
            f=open('comp_metadata/I_IOPs.txt','r')
            lines=f.readlines()
            f.close()
#            f=open('comp_metadata/I_IOPs.txt','a')
            #if line in lines:
            #    continue
            del lines
            dropsincore=0
            outerdrops=0
            #print('start filelist loop')
            # Iterate over all files.
            for filename in filelist:
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
                try:
                    r,theta=cart_to_cylindr(mlon,mlat,track,dicc['Launch Time'])
                except:
                #    print(storm,year)
                    continue

                # Continue condition only if near-inner core dropsonde.
                if r>240:
                    continue
                elif r<100:
                    dropsincore+=1
                else:
                    outerdrops+=1
            if dropsincore<13:
                if dropsincore>8 and outerdrops>5:
                    print('Missing few dropsondes')
                    print(storm,stdt,endt)
                else:
                    continue
            else:
                print('Writing to file '+storm+str(period))

            if np.isnan(get_Rmax(track,stdt,endt)):
                rmax=np.nan
            else:
                rmax=int(get_Rmax(track,stdt,endt))
            inten,Ic=periodI(storm,stdt.year,stdt,endt)

            line=storm+','+str(stdt)+','+str(endt)+','+str(inten)+','+str(Ic[0])+'\n'
            try:
                f=open('comp_metadata/I_IOPs.txt','r')

                lines=f.readlines()
                f.close()
                f=open('comp_metadata/I_IOPs.txt','a')
                if line not in lines:
                    f.write(line)
                f.close()
            except:
                f=open('comp_metadata/I_IOPs.txt','a')
                f.write(line)
                f.close()
            # Intensity-categories.
            if int(inten)>113:
                filename='comp_metadata/h4_5.txt'
            elif int(inten)>96:
                filename='comp_metadata/h3.txt'
            elif int(inten)>83:
                filename='comp_metadata/h2.txt'
            elif int(inten)>64:
                filename='comp_metadata/h1.txt'
            else:
                filename='comp_metadata/tstd.txt'
            f=open(filename,'r')
            lines=f.readlines()
            f.close()
            f=open(filename,'a')
            line=storm+','+str(stdt)+','+str(endt)+','+str(inten)+','+str(Ic[0])+','+str(dropsincore+outerdrops)+','+str(rmax)+'\n'
            if line not in lines:
                f.write(line)
            f.close()

            #IC categories

            if np.nanmean(Ic)>= 10:
                ICfile='comp_metadata/In.txt'
            elif np.nanmean(Ic)<=-10:
                ICfile='comp_metadata/We.txt'
            else:
                ICfile='comp_metadata/SS.txt'
            line=storm+','+str(stdt)+','+str(endt)+','+str(inten)+','+str(Ic[0])+','+str(dropsincore+outerdrops)+','+str(rmax)+'\n'
            f=open(ICfile,'r')
            lines=f.readlines()
            f.close()
            f=open(ICfile,'a')
            if line not in lines:
                f.write(line)
            f.close()

    # STAts COMPOSITES
#print(storm,stdt,endt,inten,Ic,dropsincore)
