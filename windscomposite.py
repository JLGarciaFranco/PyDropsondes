"""
Plot composite winds
=====================================================

This example gives the typical plot of pressure and wind speed from the best track dataset.

Specifically, the `best track <https://www.nhc.noaa.gov/data/#hurdat>`_ dataset, provides the maximum sustained winds at 10 m altitude (hereafter :math:`U_{10}`)
and the minimum surface pressure (hereafter :math:`P_{min}`).

"""
import os
import numpy as np
import scipy
import sys
import pandas as pd
import datetime
import glob
import matplotlib.pyplot as plt
from flightdata import trackandspeed
from toolbox import potential_temperature,heightavg,get_Rmax
from scipy.interpolate import griddata
ic='we/'
folder='compfiles/'+ic
filelist=glob.glob(folder+'*.txt')
filelist=np.sort(filelist)
dffilename='comp_metadata/We.txt'
df=pd.read_csv(dffilename,header=None,names=['Storm','Start','End','I',"IC","Drops","RMW"])
ri=np.arange(0,3,0.25)
Hi=np.arange(50,3200,90)
def findrmw(filename,df):
    stormperiod=filename.split('/')[-1]
    broken1=stormperiod.split('-')[0]
    broken2=stormperiod.split('-')[1]
    storm=broken1[0:-4]
    startday=int(broken1[-4:-2])
    startmin=int(broken1[-2:])

    stormdata=df[df['Storm']==storm]

    for date in stormdata['Start'].values:
        try:
            datet=datetime.datetime.strptime(date,'%Y-%m-%d %H:%M:%S')
        except:
            datet=datetime.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f')
        if datet.day == startday and datet.hour== startmin:
            stdt=datet
            endstrin=stormdata[stormdata['Start']==date]['End']
            #print(stdt)
            #print(endstrin.values)
            try:
                endt=datetime.datetime.strptime(endstrin.values[0],'%Y-%m-%d %H:%M:%S')
            except:
                endt=datetime.datetime.strptime(endstrin.values[0],'%Y-%m-%d %H:%M:%S.%f')
            break
        else:
            print(storm,date,datet.day,datet.hour,startday,startmin)

    try:
        os.system('cp /media/jlgf/Seagate\ Expansion\ Drive/FlightData/v1.1/*'+str(datet.year)+'*'+storm.upper()+'* '+'Data/'+str(datet.year)+'/'+storm+'/')
        track=trackandspeed(storm,str(datet.year))
    except:
        print(broken1,df['Storm'],filename)
        quit()
    return track,stdt,endt
#bigR=np.empty()
bigv_tang=[]
bigH=[]
for index,filename in enumerate(filelist):
    print(filename)
    track,start,endt=findrmw(filename,df)
    rmaxi=get_Rmax(track,start,endt)
    matrix=np.genfromtxt(filename)
    r=matrix[:,0]
    H=matrix[:,5]
    u_radial=matrix[:,11]
    v_tang=matrix[:,10]
    x=matrix[:,12]
    y=matrix[:,13]
    temperature=matrix[:,6]
    thetae=matrix[:,7]
    rmaxi=int(np.nanmean(r[np.where(np.nanmax(v_tang)==v_tang)]))
    if index==0:
        bigR=r/rmaxi
        bigv_tang=u_radial
        bigH=(H)
    else:
        bigR=np.append(bigR,r/rmaxi)
        bigv_tang=np.append(bigv_tang,u_radial)
        bigH=np.append(bigH,H)
    if index>7:
        break
    print(bigv_tang.shape)
    #    track=trackandspeed(storm,year)
#print(bigR)
    ri=np.arange(0,200)
    mean_radial=scipy.interpolate.griddata((r,H), u_radial,(ri[None,:], Hi[:,None]),method='linear')
    #r=r/rmaxi
    mean_azi=scipy.interpolate.griddata((r,H),v_tang, (ri[None,:], Hi[:,None]),method='linear')
    #mean_radial=scipy.interpolate.griddata((r,H),u_radial, (ri[None,:], Hi[:,None]),method='linear')
    mean_temp=scipy.interpolate.griddata((r,H),temperature, (ri[None,:], Hi[:,None]),method='linear')
    #indices = np.logical_not(np.logical_or(np.isnan(thetae), np.isnan(r)))
    mean_vert=scipy.interpolate.griddata((r,H),thetae, (ri[None,:], Hi[:,None]),method='linear')
    plt.figure(figsize=(13,13))
    plottingdictionary={"Radial wind":mean_radial,"Azimuthal wind":mean_azi,"Potential temperature":mean_temp,r"$\theta_e$":mean_vert}
    colormaps=['seismic','rainbow','coolwarm','gist_rainbow']
    labels=[r'm s$^{-1}$',r'm s$^{-1}$','K',r'K']
    spacing=[4,5,1.5,2]
    counter=0
    for variable in plottingdictionary.keys():
        field=plottingdictionary[variable]
        print(variable)
        if variable=='Radial wind':
        	maxi=-np.nanmin(field)+4
        	mini=np.nanmin(field)
        elif variable=='Azimuthal wind':
        	mini=0
        	maxi=np.nanmax(field)+5
        elif variable=='Potential temperature':
        	maxi=int(np.nanmax(field))+2
        	mini=int(np.nanmin(field))-1
        else:
        	maxi=int(np.nanmax(field))+2
        	mini=int(np.nanmin(field))
        ax=plt.subplot(221+counter)
        CS=plt.contourf(ri,Hi,field,cmap=colormaps[counter],levels=np.arange(int(mini),int(maxi),spacing[counter]))
        if 'wind' not in variable:
        	plt.xlabel('Radius [km] ',fontsize=16)
        if variable=='Radial wind' or variable=='Potential temperature':
        	plt.ylabel('Height [m] ',fontsize=16)
        #	ax.scatter(r,H,color='black',s=1)
        plt.title(variable,fontsize=18)
        plt.xlim([0,200])
        plt.ylim([0,2400])
        plt.colorbar(CS,label=labels[counter])
        counter+=1

    plt.tight_layout()
    plt.suptitle(filename[0:5],fontsize=19)
    #plt.savefig('figs/'+storm+'crossect'+str(sdate)+'.png')
    #plt.show()
    plt.close()
ri=np.arange(0,3,0.1)
mean_azi=scipy.interpolate.griddata((bigR,bigH),bigv_tang, (ri[None,:], Hi[:,None]),method='linear')
plt.contourf(ri,Hi,mean_azi,cmap='seismic')
plt.colorbar()
plt.show()
