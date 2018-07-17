"""
.. module:: main_interface

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
from processing import *
# I
# get the track from the fligh-level data.

from cape import cape
# Module for plotting track.
from plotrack import plotrack
# Module for plotting dropsonde location and drift.
#from plot_intensity_series import intensityplot
from plot_drift_dropsondes import *
# Allocate directory where files will be saved.
from flightdata import trackandspeed
# Estimate cape-measurements.

home=os.getcwd()

print(home)

figdir=home+'/figs/'

# Begin interface.
os.system('clear')
# Greeting message.
for i in range(12):
	print("")
print("          \t \t                      _     _   _   _   _  _  _        _   _  _ " )
print("          \t \t                     |_|\/ | \ |_| | | |_||_ | | |\ | | \ |_ |_                            ")
print("          \t \t                     |  /  |_/ | \ |_| |   _||_| | \| |_/ |_  _|                          ")
print("          \t \t              ============================================================                           ")
print("    _____           _                   ")
print("     | |       _   | \   _   _   _   __   _        _   _")
print("     | |  |_| |_   |  | |_| | | |_| |__  | | |\ | | \ |_      |  A visualization tool for Tropical Cyclone studies. ")
print("     | |  | | |_   |_/  | \ |_| |   ___| |_| | \| |_/ |_      |  Documentation: ")
print("        |  |    _       _         __  _  ___    _              |  email for help: ee17jlgf@leeds.ac.uk ")
print("        |  | | |_  | | |_|  |   |  / |_|  |  | | | |\ |        |  |  ")
print("         \/  |  _| |_| | |  |__ | /_ | |  |  | |_| | \|        |  |  |  Version 1.0 (2018-05-09 12:00 UTC) ")
print("         ___  _  _           ___      ")
print("          |  | || | |   |/ |  |                                  |  x86_64-linux-gnu ")
print("          |  |_||_| |__ |\ |  |     ")
print("\t -----------------------------------------------------------------------------------------------------------------------------------------------------------")
print("\t \t Welcome to the Python3.5 Interface to analyze, process and mostly visualize Tropical Cyclones using the Dropsonde Observations.")
print(" ")
print("\t \t For the Masters of Research Project, several scripts were designed to visualize very interesting features of Tropical Cyclones")
time.sleep(1)
#os.system("clear")
for i in range(12):
	print("")
print("\t \t \t \t I do not wish to get ahead of myself but you might end-up enjoying this way too much.")
print(" ")
answer=input("\t \t \t First, if you are unsure about installed packages in your distro would you like me to install the main packages used by PyDropsondes? [yes/no]  ")
print(" ")
if answer=='yes':
	print("\t \t \t \t Ok... installing ... ")
	os.system('bash install.sh')
	print(" Installation complete ---")
	time.sleep(1)
	os.system('clear')
else:
	print("\t \t \t \t Seems like I'll keep going, if it fails, look up the docs and the install.sh file.")
time.sleep(1)
print("\t \t \t \t The first piece of information I require is the storm's name and year. ")
print(" ")
time.sleep(1)
# Ask the user for the storm.
year=input("\t \t \t \t What year is your storm? (integer from 1996-2012) __ ")
print(" ")
time.sleep(1)
# Provide a list of the storms in the dataset for that year. .
print(" \t \t \t \t These are the storms that we have available")
print(" ")
time.sleep(1)
os.system('ls '+home+'/Data/'+year+'/')

print(" ")
storm=input("\t \t \t \t What storm are you looking for? __ ")

print(" ")
print("\t \t \t \t Great, you have selected Tropical Cyclone "+storm+" from the year of "+year)
time.sleep(2)
directorynames=['/gps.qc.eol/GIV/','/ublox.qc.eol/GIV/','/gps.qc.eol/P-3.43/','/ublox.qc.eol/P-3.43/','/gps.qc.eol/P-3.42/','/ublox.qc.eol/noaa.P-3/','/ublox.qc.eol/P-3.42/','/gps.qc.eol/USAF/','/ublox.qc.eol/USAF/']
#directorynames=['/gps.qc.eol/P-3.42/']
if int(year)>2012:
    downloadtype='avp'
else:
    downloadtype='radazm'
if downloadtype == 'radazm':
	filelist=[]
	for direct in directorynames:
		filelist=filelist+glob.glob(home+'/Data/'+year+'/'+storm+direct+'*')
	filelist=glob.glob(home+'/Data/'+year+'/Hawk/*')+filelist
elif downloadtype =='avp':
    filelist=glob.glob(home+'/Data/'+year+'/'+storm+'/*')
        #Define directory for output-figures (to be used by other scripts)
os.system('mkdir '+figdir+storm)
os.system('mkdir '+figdir+storm+'/planviews')
os.system('mkdir '+figdir+storm+'/axisym')
print(" ")
print('\t \t \t Retrieving track ...')

try:
	track=trackandspeed(storm,year)
except:
	print("It would seem that we do not have the track file for your storm")
#	quit()
#filelist=np.sort(filelist)
#os.system('mkdir ../figs/'+storm)

#pltleg(filelist,downloadtype,storm,track)
time.sleep(4)


os.system('clear')

print("\t \t \t We have the following activities and plots that you can do, which one would you like to perform?")
print(" ")
print("\t \t \t \t a) Storm Track")
print("\t \t \t \t b) Dropsonde drift")
print("\t \t \t \t c) Plan views and axisymmetric cross sections of cylindrical winds and temperature fields.")
print("\t \t \t \t d) Gradient wind")
print("\t \t \t \t e) Kinematic fields")
print("\t \t \t \t f) CAPE distribution")
print(" ")
print(" ")
option=input("Select a letter only __ ")
print(" ")
if option not in ['a','b','c','d','e','f','g']:
	print("I am sorry, it seems you have entered an incorrect option. Shall we try again?")
	option=input("Select a, b, c, d, e or f from the above menu.")
	if option not in ['a','b','c','d','e','f','g']:
		quit()

elif option=='a':
	print("WARNING  ----------------------------------------------")
	time.sleep(1)
	print("The option you selected requires the use of Basemap module. ")
	print("See documentation of the Plot track example on the Sphinx pages, Basemap installation and make use of the attached folder basemap-1.1.0, alternatives include the installation of Basemap through the conda package.")
	plotrack(storm,year)
elif option=='b':

	show=input("Now I only need one more piece of information, would you like me to show the results of you request on screen? (they will be saved on the figure folder in any case) [yes/no]?")
	if show =='yes':
		show=True
	else:
		show=False
	plotdrift(filelist,track,storm,show)
else:
	if len(glob.glob('outfiles/'+storm+'*'))<=0:
		print("Processing files ( this should take about 10 min.)")
		read_process_write(filelist,track,storm)
		#read_process_write(filelist,storm)

	show=input("Now I only need one more piece of information, would you like me to show the results of you request on screen? (they will be saved on the figure folder in any case) [yes/no]?")
	print("-Files processed -")
	if option == 'c':
		command='plot_winds_temp.py'
	elif option == 'd':
		command='plot_gradwindbalanc.py'
	elif option == 'e':
		command='plot_kinematic3d.py'
	elif option== 'f':
		cape(filelist,storm,track,show)
		quit()
	elif option== 'g':
		intensityplot(storm,year)
	outflist=glob.glob('outfiles/'+storm+'*')

	for filename in outflist:
		string=filename.split('/')[-1].split('.')[0]
		period=string[len(storm):]
		print('python3 '+command+' %s %s' % (storm,period))
		os.system('python3 '+command+' %s %s %s %s' % (storm,period,filename,show))
print("Routine finished -------- please find your results under figs/"+storm)
