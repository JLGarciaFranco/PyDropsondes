# -*- coding: utf-8 -*-
"""
CAPE and CIN computation
-----------------------------

The functions found below are completely random and might no be related with one another.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os
#metpy module import
from metpy.calc import saturation_mixing_ratio
from metpy.plots import add_metpy_logo, SkewT
import metpy.calc as mpcalc
from metpy.units import units,concatenate
# import all from toolbox
from toolbox import *


def cape(filelist,storm,track):
    #Sort filelist.
    filelist=np.sort(filelist)

    # Get sampling periods (this will be a dictionary). See the toolbox
    print('Retrieving sampling periods')
    sampleperiods=getsamplingperiods(filelist,3.)

    # Iterate over all sampling periods.
    for sampindex,periodskey in enumerate(sampleperiods):

        #Allocate starting (stdt) and ending date (endt). Remeber dt is the convetional short-name for date.
        stdt=periodskey
        endt=sampleperiods[periodskey]

        # Define sampling period string
        period=str(stdt.hour)+'_'+str(stdt.day)+'-'+str(endt.hour)+'_'+str(endt.day)

        # Create new-empty lists.
        lats=[]
        lons=[]
        xs=[]
        ys=[]
        capes=[]
        cins=[]

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

            nump=np.genfromtxt(filename,skip_header=16,skip_footer=0)
            temperature=clean1(nump[:,5])
            pressure=clean1(nump[:,4])
            Height=clean1(nump[:,13])
            if np.nanmax(Height)<3500:
                continue
            #Clean for cape
            RelH=clean1(nump[:,7])
            lon=clean1(nump[:,14])
            lat=clean1(nump[:,15])
            lon=clean1(lon)
            lat=clean1(lat)
            mlon=np.nanmean(lon)
            mlat=np.nanmean(lat)
            RH=RelH/100
            T,P,rh,dz=cleanforcape(temperature,pressure,RH,Height)

            #Metpy set-up
            T=np.flip(T,0)
            rh=np.flip(rh,0)
            p=np.flip(P,0)
            dz=np.flip(dz,0)
            p=p*units.hPa
            T=T*units.celsius


            mixing=rh*mpcalc.saturation_mixing_ratio(p,T)
            epsilon=0.6219800858985514
            Tv=mpcalc.virtual_temperature(T, mixing,
                                      molecular_weight_ratio=epsilon)
            dwpoint=mpcalc.dewpoint_rh(T, rh)

            blh_indx=np.where(dz<500)
            try:
                parcelprofile=mpcalc.parcel_profile(p,np.nanmean(T[blh_indx])*units.celsius,mpcalc.dewpoint_rh(np.nanmean(T[blh_indx])*units.celsius, np.nanmean(rh[blh_indx]))).to('degC')
                Tv_parcelprofile=mpcalc.virtual_temperature(parcelprofile, mixing,
                                          molecular_weight_ratio=epsilon)
                cape,cin=cape_cin(p,Tv,dwpoint,Tv_parcelprofile,dz,T)
            except:
                continue

            plotskewT=True
            if plotskewT==True:

                os.system('mkdir figs/skewt')
                fig = plt.figure(figsize=(9, 9))
                skew = SkewT(fig, rotation=45)
                skew.ax.set_ylim(1000, 100)
                skew.ax.set_xlim(-40, 60)

                skew.plot(p, dwpoint, 'g',label=r'$T_{dp}$')
                skew.plot(p, Tv, 'r',label=r'$T_v$')
                plt.text(-120,120,str(np.around(cape,2)),fontsize=14,fontweight='bold')

                # Plot the data using normal plotting functions, in this case using
                # log scaling in Y, as dictated by the typical meteorological plot
                skew.plot(p,Tv_parcelprofile,'k',label=r'$T_{v env}$')
                skew.shade_cin(p, T, parcelprofile,label='CIN')
                skew.shade_cape(p, Tv, Tv_parcelprofile,label='CAPE')
                skew.plot_dry_adiabats()
                skew.plot_moist_adiabats()

                plt.legend()
                plt.title(storm + ' on' + period,fontsize=14)
                plt.savefig('figs/skewt/'+storm+str(dicc['Launch Time'].time())+'.png')
                #plt.show()
                plt.close()

            r,theta=cart_to_cylindr(mlon,mlat,track,dicc['Launch Time'])
            if not(np.isnan(r)) and not(np.isnan(theta)) and not(np.isnan(cape.magnitude)):
                xs.append(r*np.cos(theta))
                ys.append(r*np.sin(theta))
                capes.append(cape.magnitude)
                cins.append(cin)

            fig = plt.figure(figsize=(13, 9))
            plt.scatter(xs,ys,c=np.asarray(capes),cmap='jet')
            for i,xi in enumerate(xs):
                plt.text(xi,ys[i]+10,str(np.around(capes[i],1)))

        plt.colorbar(label=r"$J/kg$')
        plt.scatter(0,0,marker='v',s=100,color='black')
        plt.grid()
        plt.xlabel('X distance [km]')
        plt.ylabel('Y distance [km]')
        plt.title('CAPE distribution for '+storm+' on '+period,fontsize=14)
        plt.savefig('figs/cape'+storm+period+'.png')
        #plt.close()
