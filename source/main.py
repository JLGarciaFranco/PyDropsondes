"""
The all powerful, yet simple interface.
----------------------------------------

This script is in charge of interacting with the user, calling other scripts and bug-control.
It will be the user's first and last impression and, as such, the script was written paying attention to both form and manner.
"""

# Import important modules.
# Numpy is the main mathematical and array-handler module in Python.
import numpy as np
from scipy.interpolate import griddata
# Module to interact with operating system.
import os
# Module to have dates and times.
import datetime
import time
#
import glob
import sys
# Import axisymmetrical plots module.
from axisym import *
# I
from hcrossections import *
# Import plan-views module.
from planviews import *
# get the track from the fligh-level data.

from cape import cape
# Module for plotting track.
from plotrack import plotrack
# Module for plotting dropsonde location and drift.

from plotleg import *
from plot_drift_dropsondes import *
# Allocate directory where files will be saved.
from flightdata import trackandspeed
# Estimate cape-measurements.

home=os.getcwd()
print(home)

figdir=home+'figs/'

# Begin interface.
os.system('clear')
# Greeting message.
for i in range(12):
	print("")
print("    _____           _                   ")
print("     | |       _   | \   _   _   _   __   _        _   _")
print("     | |  |_| |_   |  | |_| | | |_| |__  | | |\ | | \ |_      |  A visualization tool for Tropical Cyclone studies. ")
print("     | |  | | |_   |_/  | \ |_| |   ___| |_| | \| |_/ |_      |  Documentation: ")
print("        |  |    _       _         __  _  ___    _              |  email for help: ee17jlgf@leeds.ac.uk ")
print("        |  | | |_  | | |_|  |   |  / |_|  |  | | | |\ |        |  |  ")
print("         \/  |  _| |_| | |  |__ | /_ | |  |  | |_| | \|        |  |  |  Version 0.3 (2018-04-18 00:58 UTC) ")
print("         ___  _  _           ___      ")
print("          |  | || | |   |/ |  |                                  |  x86_64-linux-gnu ")
print("          |  |_||_| |__ |\ |  |     ")
print("\t -----------------------------------------------------------------------------------------------------------------------------------------------------------")
print("\t \t Welcome to the Python3.5 Interface to analyze, process and mostly visualize Tropical Cyclones using the Dropsonde Observations.")
print(" ")
print("\t \t For the Masters of Research Project, several scripts were designed to visualize very interesting features of Tropical Cyclones")
time.sleep(1)
os.system("clear")
for i in range(12):
	print("")
print("\t \t I do not wish to get ahead of myself but you might end-up enjoying this way too much.")
print("\t \t   While official names for an interface might be juvenile, you may call the Dropsonde Observations of Tropical Cyclones Visualization Tool, or just ")

# Ask the user for the storm.
year=input("What year is your storm? ")
# Provide a list of the storms in the dataset for that year. .
os.system('ls '+home+'/Data/'+year+'/')
storm=input("What storm are you looking for? ")
directorynames=['/gps.qc.eol/GIV/','/ublox.qc.eol/GIV/','/gps.qc.eol/P-3.43/','/ublox.qc.eol/P-3.43/','/gps.qc.eol/P-3.42/*radazm','/ublox.qc.eol/noaa.P-3/','/ublox.qc.eol/P-3.42/']
if int(year)>2012:
    downloadtype='avp'
else:
    downloadtype='radazm'
if downloadtype == 'radazm':
	filelist=[]
	for direct in directorynames:
		filelist=filelist+glob.glob(home+'/Data/'+year+'/'+storm+direct+'*')
	filelist=glob.glob(home+'/Data/'+year+'/'+storm+'/Hawk/*')+filelist
elif downloadtype =='avp':
    filelist=glob.glob(home+'/Data/'+year+'/'+storm+'/*')
        #Define directory for output-figures (to be used by other scripts)
os.system('mkdir '+figdir+storm)
print('Retrieving track')
track=trackandspeed(storm,year)
filelist=np.sort(filelist)
#os.system('mkdir ../figs/'+storm)

#pltleg(filelist,downloadtype,storm,track)

print("Processing file")
print(" We have the following activities and plots that you can do, which one would you like to perform?")
print(" a) Storm Track")
print(" b) Dropsonde drift")
print(" c) Plan views and axisymmetric cross sections of cylindrical winds and temperature fields.")
print(" d) Gradient wind")
print(" e) Kinematic fields")
#read_process_write(filelist,track,storm)
#plotdrift(filelist,track,storm)
cape(filelist,storm,track)
#plotuv(filelist,downloadtype,storm,track)
#plotrack(storm,year)
