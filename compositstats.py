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
filelist=['In','We','SS']
#filelist=['tstd','h1','h2','h3','h4_5']
#filelist=['I_IOPs']
plt.figure(figsize=(9,14))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
labels=['TS/TD','Cat. 1','Cat. 2','Cat. 3','Cat. 4,5','Cat. 5']
labels=['Intensifying','Weakening','Steady-state']
for counter,name in enumerate(filelist):
    filename='comp_metadata/'+name+'.txt'
    df=pd.read_csv(filename,header=None,names=['Storm','Start','End','I',"IC","Drops","RMW"])
    name=labels[counter]
    newlist=[]
    #print(df['Drops'].sum(),name)
    #continue
#    for xxx in df['IC'].values:
#        if '[' in xxx:
#            xxx=xxx.split('[')[1][0:-1]

#        if len(xxx)>5:
#            xxx=xxx[0:5]
#        print(xxx)
#        newlist.append(float(xxx))
#    df['IC']=newlist
#print(df)
#for indice in df.index:
#    series=df.loc[indice]
    #df['I'].hist(bins=[0,64,83,96,113,170])
    #df['RMW'].hist(bins=[0,20,30,40,55,75,100])
#    storm=series['Storm']
#    icstring=series['IC'][1:-2]
    #Ic=int(series['IC'])
    #print(Ic)
    rmwdict={'0-20':len(df[df['RMW']<20]),'20-30':len(df[(df['RMW']>20)&(df['RMW']<30)]),'30-40':len(df[(df['RMW']>30)&(df['RMW']<40)]),'40-50':len(df[(df['RMW']>40)&(df['RMW']<50)]),'50-65':len(df[(df['RMW']>50)&(df['RMW']<65)]),'65-80':len(df[(df['RMW']>65)&(df['RMW']<80)]),'>80':len(df[(df['RMW']>80)&(df['RMW']<110)])}
    #intensdict={'TS/TD':len(df[df['I']<64]),'Cat. 1':len(df[(df['I']>64)&(df['I']<83)]),'Cat. 2':len(df[(df['I']>83)&(df['I']<96)]),'Cat. 3':len(df[(df['I']>96)&(df['I']<113)]),'Cat. 4,5':len(df[(df['I']>113)])}
    #icdict={'We':len(df[df['IC']<=-10]),'SS':len(df[(df['IC']>-10)&(df['IC']<10)]),'In':len(df[df['IC']>=10])}
    #rmwdict=intensdict
    #rmwdict=icdict
    names = list(rmwdict.keys())
    print(names)
    values = list(rmwdict.values())
    print(values)
    #quit()
    if counter==0:
        ax=plt.subplot(311+counter)
        ax1=ax
    else:
        ax=plt.subplot(311+counter,sharex=ax1)
    ax.text(-0.5,8.5,name,bbox=props,fontweight='bold',fontsize=16)
    #plt.xlim([45,130])
    #plt.xticks([50,70,90,105,125],['TS/TD','H-1','H-2','H-3','H-4,5'])
    #plt.yticks(np.arange(0,30,5,dtype=int))
    #plt.xticks([10,25,35,50,67.5,92.5],['0-20','20-30','30-40','40-55','55-75','>75'])
    #plt.xlim([0,100])
    ax.set_yticks(np.arange(0,21,2))
    ax.bar(names,values,color='red')
    ax.set_ylim([0,10])
    ax.grid(alpha=0.5)
    if counter==1:
        ax.set_ylabel('Number of IOPs',family='italics',fontsize=19)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    if counter<2:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
    if counter==2:
        ax.set_xlabel(r'RMW [km]',fontsize=19)
#plt.suptitle('I Composites Size Histogram',fontsize=21)
    #I=np.nanmean(int(series['I']))
#plt.tight_layout()
plt.show()
