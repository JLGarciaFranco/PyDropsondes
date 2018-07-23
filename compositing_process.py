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
home=os.getcwd()
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
omega=7.2921*(10**-5)
filename='comp_metadata/I_IOPs.txt'
df=pd.read_csv(filename,header=None,names=['Storm','Start','End','I',"IC"])
#print(df)
for indice in df.index:
    series=df.loc[indice]

    storm=series['Storm']
#    icstring=series['IC'][1:-2]
    Ic=int(series['IC'])
    print(Ic)
    I=np.nanmean(int(series['I']))
    if int(I)>113:
        folder='/h4_5/'
    elif int(I)>96:
        folder='/h3/'
    elif int(I)>83:
        folder='/h2/'
    elif int(I)>64:
        folder='/h1/'
    else:
        folder='/tstd/'

    if Ic>=10:
        ifolder='/in/'
    elif Ic<=-10:
        ifolder='/we/'
    else:
        ifolder='/ss/'

    try:
        start=datetime.datetime.strptime(series['Start'],'%Y-%m-%d %H:%M:%S')
    except:
        start=datetime.datetime.strptime(series['Start'],'%Y-%m-%d %H:%M:%S.%f')
    try:
        end=datetime.datetime.strptime(series['End'],'%Y-%m-%d %H:%M:%S')
    except:
        end=datetime.datetime.strptime(series['End'],'%Y-%m-%d %H:%M:%S.%f')
    print(start,type(start),type(end))
    year=str(start.year)
    if start.hour >= 10:
        hh=str(start.hour)
    else:
        hh='0'+str(start.hour)
    if start.day >=10:
        dd=str(start.day)
    else:
        dd='0'+str(start.day)
    if end.hour<10:
        fhh='0'+str(end.hour)
    else:
        fhh=str(end.hour)
    if end.day<10:
        fdd='0'+str(end.day)
    else:
        fdd=str(end.day)
    period=dd+hh+'-'+fdd+fhh

    print(storm,period)

    outfile1='compfiles'+folder+storm+period+'.txt'
    outfile2='compfiles'+ifolder+storm+period+'.txt'
    if len(glob.glob('compfiles'+folder+storm+period+'*'))>0 and len(glob.glob('compfiles'+ifolder+storm+period+'*'))>0:
        continue
    directorynames=['/gps.qc.eol/GIV/','/ublox.qc.eol/GIV/','/gps.qc.eol/P-3.43/','/ublox.qc.eol/P-3.43/','/gps.qc.eol/P-3.42/','/ublox.qc.eol/noaa.P-3/','/ublox.qc.eol/P-3.42/','/gps.qc.eol/USAF/','/ublox.qc.eol/USAF/']
    filelist=[]

    for direct in directorynames:

    	filelist=filelist+glob.glob(home+'/Data/'+year+'/'+storm+direct+'*')

    filelist=glob.glob(home+'/Data/'+year+'/Hawk/*')+filelist
    os.system('cp /media/jlgf/Seagate\ Expansion\ Drive/FlightData/v1.1/*'+year+'*'+storm.upper()+'* '+'Data/'+year+'/'+storm+'/')
    track=trackandspeed(storm,year)
    # Create new-empty lists.
    lats=[]
    lons=[]
    for filename in filelist:



        # Select end-name of file by inspecting filename string. Notice how filename can change how file is read.
        if 'radazm' in filename.split('/')[-1] or 'eol' in filename.split('/')[-1]:
            datatype='radazm'
        else:
            datatype='avp'
        # Obtain properties of file, i.e., launch time and location into a dictionary (dicc).
        dicc=findproperties(filename,datatype)
        # Condition to see if current file is in sampling period.
        # Notice how if structure is constructed, condition finds times outside of sampling period and
        # if found outside the sampling period, continue to next file.
        if dicc['Launch Time']<start or dicc['Launch Time'] > end:
            continue

        # Allocate reading parameters.
        if datatype =='avp':
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
        elif datatype == 'radazm':
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
        r,theta=cart_to_cylindr(mlon,mlat,track,dicc['Launch Time'])

        # Continue condition only if near-inner core dropsonde.
        if r>240:
            continue
        # print filenames used.
        #print(filename)

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
        ustorm,vstorm=stormu(u[1],v[1],dicc['Launch Time'],track[3])

        # Iterate over vectores, we use longitude vector for simplicity but all vectors (height, u_speed, etc) have the same length so for
        # syntax can be used with either vector.
        for j,longi in enumerate(lon):

            # Get closest date using time arrays. Define datetime object.
            date=datetime.datetime(dicc['Launch Time'].year,dicc['Launch Time'].month,dicc['Launch Time'].day,int(hours[j]),int(minutes[j]),int(seconds[j]))

            # Check for nan, if so continue, if observation cannot be located, then all points are useless.
            if np.isnan(longi) or np.isnan(lat[j]):
                continue

            # Try and obtain radius and azimuth for current observation.
            try:
                r,theta=cart_to_cylindr(longi,lat[j],track,date)
            except:
                continue

            # Correct u_speed, make it storm relative.
            u_speed[j]=u_speed[j]-ustorm
            v_speed[j]=v_speed[j]-vstorm

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
    os.system('cp tempjulia.txt '+outfile1)
    os.system('cp tempjulia.txt '+outfile2)

    print('starting plotting sequence')
    # Call extra processing or plotting routines.
    #print('python3 3Dfields.py %s %s' % (storm,end.time()))
#            os.system('python3 gradientwind.py %s %s' % (storm,period))

    # Remove temporary file.
    os.system('rm temp_axisym.txt')
