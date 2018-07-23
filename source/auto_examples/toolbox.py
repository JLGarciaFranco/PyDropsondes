# -*- coding: utf-8 -*-
"""
Elementary tools for the DOT
-----------------------------

The functions found below are completely random and might no be related with one another.

"""
import numpy as np
import datetime
import pandas as pd
from scipy.interpolate import griddata
from scipy import interpolate
from math import sin, cos, sqrt, atan2, radians
import metpy.calc as mpcalc
from metpy.calc import find_intersections
from metpy.units import units,concatenate
import matplotlib.pyplot as plt
def distance(lat1,lon1,lat2,lon2):
	r"""Estimate the distance in kilometers between two points.

	Given the latitude and longitude of two points in Earth (x1,y1), (x2,y2)
	This function computes the distance between those two points using formulaes of spheric geometry.


	*Parameters*

	lat1 : `float`
	    Latitude of starting point [radians]
	lon1 : `float`
	    Longitude of starting point [radians]
	lat2 : `float`
	    Latitude of final point [radians]
	lon2 : `float`
	    Longitude of final point [radians]

	*Returns*

	r : `float`
	     The geometric distance between two points on the surface of a sphere

	The radius *r* is given using, first, the
	`Haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_ and it was written to match the script from `Rosetta <https://rosettacode.org/wiki/Haversine_formula#Python>`_ to great extent.

	.. math:: a = sin^2\bigg(\frac{\Delta \varphi}{2}\bigg) +cos \varphi _1 \cdot cos \varphi _2 \cdot sin^2\bigg(\frac{\Delta \lambda}{2}\bigg)

	determined by:

	* :math:`a` Great-circle distance between to points.
	* :math:`\Delta \varphi=\varphi _2-\varphi _1` Latitude difference between two points.
	* :math:`\varphi_1` Latitude of starting point [radians].
	* :math:`\varphi_2` Latitude of final point [radians].
	* :math:`\Delta \lambda=\lambda _2-\lambda _1` Longitude difference between two points.

	To convert the great circle distance given by the Haversine formula in a sphere of radius :math:`R` (radius of the Earth) to a distance :math:`r`:

	.. math:: r=R\cdot 2 arctan\bigg(\frac{\sqrt{a}}{\sqrt{1-a}}\bigg)

	Using the previous mathematical treatment, you could use this functions as desribed below:


	*Examples*

	>>> from toolbox import distance
	>>> from math import radians
	>>> print(distance(radians(45),radians(80),radians(59),radians(82)))
	1563.051273807386


	"""
	# approximate radius of earth in km
	R = 6373

	# delta lambda y varphi
	dlon = lon2 - lon1
	dlat = lat2 - lat1

	# greater circle distance
	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2

	# split operation to get radius
	factor= 2 * atan2(sqrt(a), sqrt(1 - a))

	# distance is the radius of the Earth times the parameter c
	distance = R * factor

	# Return output.
	return distance
# Get track function of storm and year.
def getrack(storm,year):
    r"""
    Obtain the track of a storm specifying the year since storm names can be repeated. .

    This function obtains and returns all track parameters, location, intensity and relative speed.

    *Parameters*

    storm : `string`
        Name of storm.
    year : `string`
        Year of storm.

    *Returns*

    longitudes :`np.ndarray (type=np.float)`
    latitudes:`np.ndarray(type=np.float)`
    windspeeds:`np.ndarray(type=np.float)`
    Intensity_label: `list`
        List of Intensity Labels, e.g., (H) for Hurricane.
    df.index:`pandas.DataFrame.Index`
        List of dates that correspond to the dates and times of the previous arrays.

    *Examples*

    >>> from toolbox import getrack
    >>> track=getrack('Isabel','2003')


    """
def clean1(vec):
	r"""
	Clean values from arrays that are NANs but written as any of the following:

	*9999.0,999.0,99.0,99999.0,-999.,-9999.,-99. *

	This function simply replaces this values for Not a Number (NAN) values that Python can recognize.

	*Parameters*

	vec : `np.array`
	    Array to clean

	*Returns*

	vec :`np.ndarray (type=np.float)`

	*Examples*

	>>> from toolbox import clean1
	>>> newpressure=clean1(pressure)


	"""
	falseval=[9999.0,999.0,99.0,99999.0,-999.,-9999.,-99.]
	for i,v in enumerate(vec):
		if v in falseval:
			vec[i]=np.nan
	return vec
def cleanp(vec):
	for i,v in enumerate(vec):
		if v < 500 or v>1100:
			vec[i]=np.nan
	return vec
def clean2(vec):
	for i,v in enumerate(vec):
		if np.isnan(v) and i < len(vec)-1:
			if not(np.isnan(vec[i-1])) and not(np.isnan(vec[i+1])):
				vec[i]=(vec[i-1]+vec[i+1])/2.0
			elif not(np.isnan(vec[i-1])) and i>4:
				vec[i]=np.nanmean(vec[i-4:i])
	return vec
def cleanu(vec):
	for i,v in enumerate(vec):
		if np.abs(v) > np.abs(np.mean(vec)+3*np.std(vec)):
			vec[i]=np.nan
		if np.abs(v-vec[i-1]) > np.abs(1.1*np.std(vec)):
			vec[i]=np.nan
	dif=np.diff(vec)
	ii=np.where(np.abs(dif)>1)
	for i in ii:
		vec[i+1]=np.nan
	return vec
def dp_dr(pressure, radius):
	r"""Calculate the pressure gradient on the radius coordinate.

	This function estimates the gradient of the function :math:`p` along the radial coordinate, *i.e.*:

	.. math:: \frac{\partial p}{\partial r}\sim \frac{\Delta p}{\Delta r}=\frac{P_2-P_1}{r_2-r_1}

	*Parameters*

	pressure : `np.ndarray (type=np.float)`
	    The total interpolated atmospheric pressure in milibars.
	radius : `np.ndarray (type=np.float)`
	    Radius distance to the center of the storm in [km].

	*Returns*

	`np.ndarray (type=np.float)`


	*Examples*

	>>> from toolbox import dp_dr
	>>> dpdr=dp_dr(pressure,r)))

	"""
	#Convert radius to [m]
	radius=radius*1000.
	# Convert pressure to Pascals from hPa
	pressure=pressure*100

	# Empty lists to fill: derivate and new radius vector.
	dpdr=[]
	newr=[]


		#  get first_derivative loop
	for index,p1 in enumerate(pressure):

		# Estimate pressure difference, p2-p1
		deltap=pressure[index+1]-p1
		# r1=radius[index] and r2[index+1]
		deltar=radius[index+1]-radius[index]
		# dp_dr according to the equation above
		dpdr.append(deltap/deltar)

		# Get mean radius where derivative was estimated
		newr.append((radius[index+1]+radius[index])/2.)

		# Break condition of final value
		if index+1==len(pressure)-1:
			break
	# Return array and new radius vector
	return np.array(dpdr),np.asarray(newr)
def getgradwind(presgrad, radius,coriolis):
	r"""Calculate the gradient wind in a particular point.

	Gradient wind balance is given by

	.. math:: \frac{1}{\rho_0}\frac{\partial p}{\partial r}=\frac{V_g^2}{r}+fV_g

	where :math:`\rho_0` is a constant density, :math:`\frac{\partial p}{\partial r}` is the
	radial gradient of the pressure field, :math:`V_g` is the gradient wind and :math:`f` is the coriolis parameter.

	which is a polynomial of degree *n=2* and needs to be solved using
	`np.roots <https://docs.scipy.org/doc/numpy/reference/generated/numpy.roots.html>`_.

	*Parameters*

	presgrad : `np.float (array)`
	    Radial gradient of pressure. See also :ref:`dp_dr`
	radius : `np.float (array)`
	    Radius distance to the center of the storm in [km].


	*Returns*

	`np.ndarray (type=np.float)`
		Gradient wind vector


	*Examples*

	>>> from toolbox import dp_dr
	>>> dpdr=getgradwind(dpdr,r0,1*(10**(-4))))

	"""
	# Allocate variables
#	Density approximation
	density=1.11164
# coriolis
	f=coriolis

	# empty list to fill with vg values
	gradwind=[]

	# np. roots requires a,b,c coefficients, see docs for np.roots.
	# ax^2+bx+c=0

	for index,r0 in enumerate(radius):
		# c
		c=-(1/density)*(presgrad[index])
		# b
		b=f[index]
		# a
		a=1/r0

		# gradientwind
		r=np.roots([a,b,c])
		vg=np.nanmax(r)
		gradwind.append(vg)

	# Return array and new radius vector
	return np.array(gradwind)
def potential_temperature(pressure, temperature):
	r"""Calculate the potential temperature.

	This function originated from the Metpy Module. It was modified by JLGF for the dropsondes.
	Uses the Poisson equation to calculation the potential temperature
	given `pressure` and `temperature`.

	*Parameters*

	pressure : `pint.Quantity`
	    The total atmospheric pressure in milibars.
	temperature : `pint.Quantity`
	    The temperature in Kelvin.

	*Returns*

	`pint.Quantity`
	    The potential temperature corresponding to the temperature and
	    pressure.

	Formula:

	.. math:: \Theta = T (P_0 / P)^\kappa

	*Examples*

	>>> from metpy.units import units
	>>> from toolbox import potential_temperature
	>>> print(potential_temperature(800. * units.mbar, 273. * units.kelvin))
	<Quantity(290.96653180346203, 'kelvin')>

	"""
	#Reference pressure
	P0=1000*units.mbar
	#  specific heat at constant pressure for dry air, in J / kg / K
	cp = 1004.
	#  gas constant for dry air, in J / kg / K
	Rd = 287.
	# Kappa Rd /Cp
	kappa = Rd / cp
	return temperature * (P0 / pressure).to('dimensionless')**kappa
def findproperties(filename,database):
	"""
	Function to read-in and find the properties of a file.

	*Parameters*

	filename: `string`
	    Name of file (full path).
	database : `dictionary`
	    The temperature in Kelvin.

	**Properties**
		1. Date/Time of sounding.
		2. Location of initial sounding.
		3. Name of sounding.

	**Returns**

	diccionario:`dictionary`

	Notice the variable-name in spanish. Given the similarity of the word with the english version and the fact that
	`dictionary <https://docs.python.org/2/tutorial/datastructures.html>`_ is a reserved word in python.

	*See also*

	:meth:`getrack`

	"""
	# Determine type of file
	if database=='avp':
		# Determine where our variables of interest are in this particular file.
		indexes=[-17,-6,-15]
		# Determine datetime format according to type of sounding.
		formato='%Y/%m/%d, %H:%M:%S.%f\n'
	# Same for radazm type of file.
	elif database=='radazm':
		indexes=[2,4,5]
		formato='%Y, %m, %d, %H:%M:%S '

	# Open -read and close file to save memory.
	f=open(filename,'r')
	lineas=f.readlines()
	f.close()

	# Create dictionary and define their keys.
	diccionario={'Sounding name':' ','lon,lat,alt':' ','Launch Time':' '}

	# Select lines (l) of the name of sounding (lname), location (location) and time (ltime).
	lname=lineas[indexes[0]].split(':')
	location=lineas[indexes[1]].split(':')
	ltime=lineas[indexes[2]].split('):')
	# Possible print for user
	#print(lname,location,ltime)

	# Allocate sounding name in dictionary.
	diccionario['Sounding name']=lname[-1]

	#Split line of location to get only relevant info.
	location=location[1].split(',')
	# Allocate location of drop in dictionary.
	diccionario['lon,lat,alt']=location[1:]

	# Allocate launch time in dictionary.
	diccionario['Launch Time']=ltime[1]

	# Adjust launch time to get only the string of launch time. Initially, diccionario['Launch Time'] has a lot of white-space.
	clear_white=diccionario['Launch Time']
	# while loop to eliminate white=space ' '
	while clear_white[0]==' ':
		clear_white=clear_white[1:]
	# Change Launch time to datetime object.
	diccionario['Launch Time']=datetime.datetime.strptime(clear_white,formato)


	# Similar routine to clean sounding Name.
	clear_white=diccionario['Sounding name']

	while clear_white[0]==' ':
		clear_white=clear_white[1:]
	diccionario['Sounding name']=clear_white

	# Return styled-dictionary.
	return diccionario

def timeconversion(hhmmss):
	"""
	Function that converts a time string to extract hours minutes and seconds from a string.
	Specific for AVAPS dropsondes.

	*Parameters*

	hhmmss: `string`
	    String containing hours (hh), minutes (mm) and seconds (ss)

	**Properties**
		1. Date/Time of sounding.
		2. Location of initial sounding.
		3. Name of sounding.

	**Returns**

	hours:`int`
	seconds:`int`
	minutes:`int`

	The rationale behind this script is that most scripts work with dates and times oriented to datetime objects.
	In this case, this particular string poses difficulties to process and as such, we extract the integers of the time.

	*See also*

	:meth:`findproperties`

	"""
	# Create empty numpy arrays to be filled in processing loop.
	hours=np.zeros(len(hhmmss))
	minutes=np.zeros(len(hhmmss))
	seconds=np.zeros(len(hhmmss))

	#Processing loop iterating over all values in hhmmss in an enumerated way.
	# Index is an integer index=0,...n. and string is the value in the array.
	for index,string in enumerate(hhmmss):
		#Obtainining first value of split string.
		string=str(string).split('.')[0]
		#Condition to see if hour is less than 10, then add a zero to read in a universal format.
		# Condition is based on length of the string, for instance 12545 corresponds to hour 1, minute 25 and 45 seconds,
		# whereas 123432 has length 6, and hour is 12.
		while len(string)<=5:
			string='0'+string

		# Allocate values in string to hours, minutes and seconds.
		hours[index]=int(string[0:2])
		minutes[index]=int(string[2:4])
		seconds[index]=int(string[4:6])

	# Return tuple (3 values in one python Object, must consider when reading output from this function)
	return hours,minutes,seconds
def getleg(lats,lons,r0):
	newlats=lats[:]
	newlongs=lons[:]
	counter=0
	for i in range(len(lats)):
		#print(newlats[i-counter],newlongs[i-counter])
		del newlats[i-counter]
		del newlongs[i-counter]
	#	print(i,counter)
		r=np.corrcoef(newlongs,newlats)
		r=r[0,1]
		print(r0,r)
		#print(newlats,lats)
		if np.abs(r0)>0.935 or len(newlats)<4:
			break
		if np.abs(r) > np.abs(r0)+0.01:
			r0=r
			counter+=1
			continue
		else:
			newlats.insert(i-counter,lats[i])
			#print(i,len(lons),len(lats))
			newlongs.insert(i-counter,lons[i])
	#print(r,len(newlats))
	return newlats,newlongs,r
def getsamplingperiods(flist,tspan):
	"""
	*Sampling Periods*


	**Parameters**

	flist: `list`
	    List of path/filename objects with all times where a dropsonde was measured.
	tspan : `float`
	    Time threshold. Maximum time allowed of difference between dropsonde measurements in one period.

	This function gets the sampling periods of a storm.
	A sampling period is a continous period of dropsonde measurements that provides un-interrupted measurements
	with a maximum threshold (tspan) specified between measurements.
	For instance,

	.. code-block:: python

		getsamplingperiods(dates,3.)

	will allow 3 hours as the maximum time between measurements to define a sampling period. If, there are more than 3 hours
	between a sequence of measurements and the next, the two sequences will account for two different-sets of measurements, i.e., two different sampling periods.

	*Returns*

	sampleperiods:`dictionary`

	See also `dictionary <https://docs.python.org/2/tutorial/datastructures.html>`_ syntax.


	"""
	# Empty list object.
	dates=[]
	# Iteration over files.
	for filename in flist:
		# Establish type of file.
		if 'radazm' in filename.split('/')[-1] or 'eol' in filename.split('/')[-1]:
			ftype='radazm'
		else:
			ftype='avp'
		# Called property finding function.
		dicc=findproperties(filename,ftype)

		# Append to dates list.
		dates.append(dicc['Launch Time'])

	# Sort dates, to start with earliest measurment.
	dates.sort()

	# Create dictionary object.
	sampleperiods={}

	# Counter
	counter=0

	# while loop to avoid Error of large indexes to date=array
	while counter<len(dates)-2:
		# Select initial date.
		dt0=dates[counter]
		# Compute difference between this date and following date.
		dif=dates[counter+1]-dt0
		# Estimate hourly difference.
		hours, remainder = divmod(dif.seconds, 3600)

		# Second while loop
		while hours <=tspan and dif.days==0:
			# Add to counter for interation continues.
			counter+=1
			# Compute difference, notice explicit use of dates[counter] since dt0 must be saved to be dictionary key.
			dif=dates[counter+1]-dates[counter]
			# Get hours again.
			hours, remainder = divmod(dif.seconds, 3600)

			# Break sequence. Could this be commented out given the outer while-loop syntax?
			if counter+2==len(dates):
				break
			# no.

		#Add key and content to dictionary using the first and last datetime objects in this period.
		sampleperiods[dt0]=dates[counter]
		# Add to counter so next while loop
		counter+=1
	return sampleperiods
def cleanforcape(T,p,dwpoint,dz):
	"""
	This function was created to clean fields relevant for CAPE calculations from NAN values.
	In simple words, this function makes sure that all four input arguments are returned where ALL arrays have non-nan values at all points.
	This cleaning function will allow further calculations to run smoothly.

	*Parameters*

	p: `pint.Quantity`
	    The total atmospheric pressure in milibars.
	T : `pint.Quantity`
	    The temperature in Kelvin.
	dwpoint :`np.array`
		Dewpoint or relative humidity.
	dz: `np.array`
		Height array.

	To be honest, this functon only nan-checks the lists and selects all indexes where no nans are found across the arrays.
	"""
	newT=[]
	newdz=[]
	newp=[]
	newdw=[]
	for i,temp in enumerate(T):
		if np.isnan(temp) or np.isnan(p[i]) or np.isnan(dwpoint[i]) or np.isnan(dz[i]):
			continue
		else:
			print(temp,p[i],dwpoint[i])
			newT.append(temp)
			newp.append(p[i])
			newdw.append(dwpoint[i])
			newdz.append(dz[i])

	return newT,newp,newdw,newdz
def interp(H,T,tipo):
#	print(H.shape[1])
#	print(H.shape)
	if tipo == 'height':
		jump=10
	elif tipo == 'pressure':
		jump=2
	minh=0
	maxh=3000
	#print(H)
	for i in range(0,H.shape[1]):
		h=H[:,i]
		slvec=T[i,:]

		if np.nanmin(h)>minh:
			#minh=np.nanmin(h)
			slvec=np.insert(slvec,0,np.nan)
			h=np.insert(h,0,minh)
		if np.nanmax(h)<maxh:
			#maxh=np.nanmax(h)
			slvec=np.insert(slvec,-1,np.nan)
			h=np.insert(h,-1,maxh)
	#	if len(h)>H.shape[0]:
			#H=np.reshape(H,(len(h),H.shape[1]))
		#H[:,i]=h
		#T[i,:]=slvec
	#print(minh,maxh)
	hnew=np.arange(0,maxh,jump)
	tnew=np.empty([H.shape[1],len(hnew)])
	for i in range(0,H.shape[1]):
		t=T[i,:]
		h=H[:,i]
	#	print(len(h),len(t))
		#griddata((xs, ys), u, (xaxis[None,:], yaxis[:,None]), method='cubic')
		f = interpolate.interp1d(h, t,fill_value=np.nan,bounds_error=False)
		ts=f(hnew)
		tnew[i,:]=ts
		#plt.plot(T[i,:],H[:,i])
		#plt.show()
	#	plt.plot(ts,hnew,label=str(i))
	#plt.legend()
#	plt.show()
	return hnew,tnew
#def filling2(h,vec):

def refill(h,maxr):
	if len(h)<maxr:
		counti=len(h)
		while counti<maxr:
			h=np.insert(h,counti,np.nan)
			counti+=1
	return h
#Function to estimate and u and v relative to the storm montion
def plotdrift(filelist,track,storm):
	"""
	:ref:`sphx_glr_auto_examples_plot_drift.py`

	*Parameters*

	filelist: `list`
		List of files to be used.
	track: `dictionary`
		Dictionary of track information. See :meth:`getrack`
	storm: `string`
		Storm name.

	*Returns*
		Figure Object


	"""
def stormu(u,v,date,dicc):
	r"""Calculate the corresponding storm speed.

	This function reads-in the storm speed at a particular datetime.

	*Parameters*

	u : `np.float`
	    X-component of wind speed.
	v : `np.float`
	    Y-component of wind speed
	date : `datetime.datetime`
	 	Date where speed is needed.
	dicc : `dictionary`
		Track dictionary with, track, speed and corresponding datetimes.

	*Returns*

	`tuple` : (newu,newv)
	    Either returns new storm-relative x and y-component wind speeds or it can return

	Formula:

	.. math:: u_{sr}=u_{obs}-u_{storm}

	where :math:`u_{sr}` is the storm relative winds, :math:`u_{obs}` is the observed wind and :math:`u_{storm}` is the estimated
	storm wind speed.

	*Examples*

	>>> from toolbox import stormu
	>>> newu,newv=stormu(u,v,date,trackdicc)

		"""
	# Extract date array from dictionary
	dates=dicc['Datetime']

	# Convert dates to object datetime
	for i,dt in enumerate(dates):
		dates[i]=pd.to_datetime(dt)


	# Extract u and v speeds from dictionary
	us=dicc['U']
	date=date+datetime.timedelta(minutes=60)
	vs=dicc['V']

	# Empty lists for storm speeds.
	uu=[]
	vv=[]

	# Loop to find speeds in and around 1 hour of our current date (date).
	for i,dt in enumerate(dates):
		# 2 h threshold
		if dt>date-datetime.timedelta(hours=2) and dt<date+datetime.timedelta(hours=2):
			uu.append(us[i])
			vv.append(vs[i])

	# Follow formula using mean of storm speeds around 1 hour
	newu=u-np.nanmean(uu)
	newv=v-np.nanmean(vv)
	# Return new u and new v.
	return np.nanmean(uu),np.nanmean(vv)
def backtoxy(rs,thetas,u,v,trackdata):
	xs=[]
	ys=[]

	for i,r in enumerate(rs):
		xs.append(r*cos(thetas[i]))
		ys.append(r*sin(thetas[i]))
	xaxis=np.linspace(np.min(xs)-1,np.max(xs)+1,0)
	yaxis=np.linspace(np.min(ys)-1,np.max(ys)+1,100)
	uinterp = griddata((xs, ys), u, (xaxis[None,:], yaxis[:,None]), method='linear')
	vinterp = griddata((xs, ys), v, (xaxis[None,:], yaxis[:,None]), method='linear')

	return xaxis,yaxis,uinterp,vinterp
def getcenter(date,track):
	r""" Obtain longitude and latitude of storm centre.

	This function is part of a basic part of most routines as it reads-in a date and the track-dictionary and returns the best approximation
	to the stomr's centre in degrees.

	*Parameters*

	date : `datetime`
	    Date to seek storm centre.
	track : `dictionary`
	    Typical dictionary obtained from :meth:`flightdata.trackandspeed`

	*Returns*

	`tuple` : (latf,lonf)
	    Return of two floats as a tuple, the final latitude (latf) and longitude (lonf).


	This script uses two main approaches, first seeking all track data in a 20 min window to obtain an average and,
	if no track values are reported in the window, the closest track value is used.

	*Examples*

	>>> from toolbox import getcenter
	>>> centrelat,centrelon=getcenter(datetime.datetime(2003,9,14,2,0,0),trackdict)

	.. note::

		Notice this script makes use of conventional and particular syntax. For instance, timedeltas are conventionally named (td) in datetime modules.


	"""
	#Unpack data from track big dictionary into a datetime list (trackdates), latitudes (traclat), longitude (traclon) and allocates the speed dictionary (not used here).
	trackdates,traclat,traclon,dicc=track

	#Grab first date in track record.
	goodate=trackdates[0]

	# Account for date difference (track is 1 hour ahead)
	date=date+datetime.timedelta(minutes=60)

	# Select how difference or timedelta (td) in datetimes is computed, the oldest date minus the earliest date so difference is positive.
	if goodate>date:
		td=goodate-date
	else:
		td=date-goodate

	# Create empty lists to be used as arrays to get average storm centre from all close times.
	lats=[]
	lons=[]

	# Iterate over all track dates tdat.
	for index,tdat in enumerate(trackdates):
		# Condition to determine how timedelta (td) is computed.
		if date > tdat:
			newtd=date-tdat
		else:
			newtd=tdat-date

		# Timedelta threshold of 10 minutes. Notice that from the timedelta definition, this is a 20 min. window.
		if newtd<datetime.timedelta(minutes=10):
			# Append latitude and longitude in this time-step, if condition is true.
			lats.append(traclat[index])
			lons.append(traclon[index])
			# Possible user print, to see exactly how track location and times are sliced.
			#print(date,goodate,dates,traclat[i],traclon[i])

		# Condition to find closest track time to the observation time.
		# In other words, finding the minimum timedelta.
		if newtd<td:
			# Goodate is not used but it is the closest time and it is good to keep it as a reference and might be printed above.
			goodate=tdat

			# Allocate closest latitude and longitude. Variables will be overwritten when a closer time is found.
			lat=traclat[index]
			lon=traclon[index]

	# If more than 1 track time is found in the 20 min window, the track centre is defined as the mean from the list of
	# all near centre observations
	if len(lats)>=2:
		# possible user print to see the lats list.
		#print(lats)

		# Define final latitude (latf) and longitude (lonf).
		latf=np.nanmean(np.array(lats))
		lonf=np.nanmean(np.array(lons))

	# If not enough close times exist in the window, track centre is reported as the closest time.
	else:
		latf=lat
		lonf=lon

	# Return tuple of storm centre latitude and longitude .
	return latf,lonf
def cart_to_cylindr(lon,lat,track,dates):
	r"""

	* Cartesian to Cylindrical *

	Convert coordinates from a dropsonde observation in cartesian coordinates :math:`(x,y)` to a cylindrical system :math:`(r,\theta)`.

	This function is part of the main part of most scripts since Tropical Cyclones are usually depicted in cylindrical systems, not in cartesian.
	As such, extreme care was taken into this function and reader is advised to look carefully at all steps of this function.

	*Parameters*

	lon : `float`
	    Longitude of observation point in degrees.
	lat : `float`
	    Latitude of observation point in degrees.
	track: `dictionary`
		Track dictionary from :meth:`flightdata.trackandspeed`.
	dates: `datetime`
		Date and time of observation.

	*Returns*

	`tuple` : (r,theta)
	    Return of two floats as a tuple, the final radius (r) and azimuth (theta).

	While this script could be self-contained it mostly depends on two other functions :meth:`distance` and
	:meth:`getcenter`. The latter retrieves the closest storm centre and the former computes distances between two lat-lon points.

	While :meth:`distance` documents how a distance between two points on a sphere is computed :math:`(r)``, the other part of a cylindrical
	coordinate system is the azimuth :math:`\theta`.

	That computation obeys a different mathematical approach. Specifically, to estimate the `azimuth <https://en.wikipedia.org/wiki/Azimuth>`_ between two points on the surface of a sphere
	nautical term of *bearing* is used.

	The bearing or azimuth is given by:

	.. math:: \theta=\frac{\pi}{2}+arctan\bigg(sin(\Delta \lambda)cos(\varphi_1),cos(\varphi_1)sin(\varphi_2)-sin(\varphi_1)cos(\varphi_2)*cos(\Delta \lambda)\bigg)

	determined by:

	* :math:`theta`: Azimuth/bearing in a cylindrical system.
	* :math:`\Delta \varphi=\varphi _2-\varphi _1` Latitude difference between centre and observation.
	* :math:`\varphi_1` Latitude of centre [radians].
	* :math:`\varphi_2` Latitude of observation [radians].
	* :math:`\Delta \lambda=\lambda _2-\lambda _1` Longitude difference between centre and observation.


	Notice how the formula has a correcting factor which accounts for the fact that a bearing usually is measured as the angle between the north (remember this was initially computed for sailing) and the location of the destination.
	For this reason, we add :math:`=\frac{\pi}{2}` to establish a proper cylindrical system where the 0\degree is aligned with the :math:`x>0` and :math:`y=0` line.

	*Examples*

	>>> from toolbox import getcenter
	>>> centrelat,centrelon=getcenter(datetime.datetime(2003,9,14,2,0,0),trackdict)

	.. note::

		For further information, good revisions on the derviation of the `Havesine formula <https://www.math.ksu.edu/~dbski/writings/haversine.pdf>`_ and the `bearing <https://www.movable-type.co.uk/scripts/latlong.html>`_ are attached.


	"""
	# Conversion of observation points to radians.
	lat2 = radians(np.abs(lat))
	lon2 = radians(np.abs(lon))

	# Obtain enter latitudes and longitudes.
	clat,clon=getcenter(dates,track)

	# Conversion to radians of centre coordinates.
	lat1=radians(np.abs(clat))
	lon1=radians(np.abs(clon))

	# Compute distance between centre and observation point.
	r=distance(lat1,lon1,lat2,lon2)

	# Compute deltas of longitude and latitudes.
	dlon = lon2 - lon1
	dlat = lat2 - lat1

	# Compute bearing in one step.
	theta2=(np.pi/2.)+np.arctan2(np.sin(dlon)*np.cos(lat2),np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(dlon))

	# Rearrange negative bearing to be positive, theta is then a subspace of azimuths from 0 ... 2 pi.
	if theta2<0:
		theta2=2*np.pi+theta2

	# Return cylindrical coordinate system.
	return r,theta2
def reshape(lat,lon,dicc):
	lon=np.array(lon)
	lat=np.array(lat)
	newdicc={}
	p=lon.argsort()
	newlon=lon[p]
	newlat=lat[p]
	for i,longi in enumerate(newlon):
		key=str(longi)+','+str(newlat[i])
		nkey=longi
		newdicc[nkey]=dicc[key]
	distances=[0]
#	print(newlon,lon)
	lon0=newlon[0]
	lat0=newlat[0]
	ii=1
	while ii <= len(newlon)-1:
		lat2=newlat[ii]
		lon2=newlon[ii]
		r=distance(lat0,lon0,lat2,lon2)
		ii+=1
		distances.append(r)
	for i,r0 in enumerate(distances):
		if r0>400:
		#	print(r0)
			del distances[i]
			del newdicc[lon[i]]
			newlat=np.delete(newlat,i)
			newlon=np.delete(newlon,i)
#	plt.plot(distances)
	#plt.show()
	#print(distances)
	return newlat,newlon,newdicc,distances
def findvalues(z,level):
        i=0
        zi=z[i]
        index=[]
        while np.abs(zi-level)>10 or np.isnan(zi):
                i+=1
                zi=z[i]
                if i==len(z)-3:
                        return
        index.append(i)
        i+=1
        zi=z[i]

        while np.abs(zi-level)<10 or np.isnan(zi):
                index.append(i)
                i+=1
                zi=z[i]
                if i==len(z)-3:
                        return
        return index
def divergence(u, v, dx, dy):
    r"""Calculate the horizontal divergence of the horizontal wind.

    **Parameters**

    u : (M, N) ndarray
        x component of the wind
    v : (M, N) ndarray
        y component of the wind
    dx : float
        The grid spacing in the x-direction
    dy : float
        The grid spacing in the y-direction

    **Returns**

    (M, N) ndarray
        The horizontal divergence

    See Also

    :ref:`stormu`

    """
    dudx = first_derivative(u, delta=dx, axis=1)
    dvdy = first_derivative(v, delta=dy, axis=0)
    return dudx + dvdy
def equivalent_potential_temperature(pressure, temperature, dewpoint):
	r"""Calculate equivalent potential temperature.

	This calculation must be given an air parcel's pressure, temperature, and dewpoint.
	The implementation uses the formula outlined in [Bolton1980]_:

	First, the LCL temperature is calculated:

	.. math:: T_{L}=\frac{1}{\frac{1}{T_{D}-56}+\frac{ln(T_{K}/T_{D})}{800}}+56

	Which is then used to calculate the potential temperature at the LCL:

	.. math:: \theta_{DL}=T_{K}\left(\frac{1000}{p-e}\right)^k
	          \left(\frac{T_{K}}{T_{L}}\right)^{.28r}

	Both of these are used to calculate the final equivalent potential temperature:

	.. math:: \theta_{E}=\theta_{DL}\exp\left[\left(\frac{3036.}{T_{L}}
	                                          -1.78\right)*r(1+.448r)\right]

	**Parameters**

	pressure: `pint.Quantity`
	    Total atmospheric pressure
	temperature: `pint.Quantity`
	    Temperature of parcel
	dewpoint: `pint.Quantity`
	    Dewpoint of parcel

	**Returns**

	`pint.Quantity`
	    The equivalent potential temperature of the parcel

	Notes
	-----
	[Bolton1980]_ formula for Theta-e is used, since according to
	[DaviesJones2009]_ it is the most accurate non-iterative formulation
	available.

	"""
	t = temperature.to('kelvin').magnitude
	td = dewpoint.to('kelvin').magnitude
	p = pressure.to('hPa').magnitude
	e = mpcalc.saturation_vapor_pressure(dewpoint).to('hPa').magnitude
	r = mpcalc.saturation_mixing_ratio(pressure, dewpoint).magnitude
	kappa=0.2854
	t_l = 56 + 1. / (1. / (td - 56) + np.log(t / td) / 800.)
	th_l = t * (1000 / (p - e)) ** kappa * (t / t_l) ** (0.28 * r)
	th_e = th_l * np.exp((3036. / t_l - 1.78) * r * (1 + 0.448 * r))

	return th_e * units.kelvin
def _gradient(f, *args, **kwargs):
    """Wrap :func:`numpy.gradient` to handle units."""
    if len(args) < f.ndim:
        args = list(args)
        args.extend([units.Quantity(1., 'dimensionless')] * (f.ndim - len(args)))
    grad = np.gradient(f, *(a.magnitude for a in args), **kwargs)
    if f.ndim == 1:
        return units.Quantity(grad, f.units / args[0].units)
    return [units.Quantity(g, f.units / dx.units) for dx, g in zip(args, grad)]


def _get_gradients(u, v, dx, dy):
    """Return derivatives for components to simplify calculating convergence and vorticity."""
    dudy, dudx = _gradient(u, dy, dx)
    dvdy, dvdx = _gradient(v, dy, dx)
    return dudx, dudy, dvdx, dvdy

def convergence_vorticity(u, v, xvec,yvec, dim_order='xy'):
	r"""Calculate the horizontal divergence and vertical vorticity of the horizontal wind.

	Parameters
	----------
	u : (M, N) ndarray
		x component of the wind
	v : (M, N) ndarray
		y component of the wind
	dx : float
		The grid spacing or vector in the x-direction
	dy : float
		The grid spacing or vector in the y-direction

	Returns
	-------
	divergence, vorticity : tuple of (M, N) ndarrays
	The horizontal divergence and vertical vorticity, respectively

	"""
	Xgrid,Ygrid=np.meshgrid(xvec,yvec)
	print(Xgrid.shape)
	indicex=np.where((xvec<=0))[0]
	indicy=np.where(yvec<=0)[0]
	gradu = np.gradient(u, xvec*1000,yvec*1000)
	gradv =np.gradient(v,xvec*1000,yvec*1000)
	dudy,dudx=gradu
	dvdy,dvdx=gradv

	return (dudx + dvdy), (dvdx - dudy)
def reassemble(r,matrix,H):
	newr=np.sort(r)
	newmatrix=np.zeros(matrix.shape)
	for i,r0 in enumerate(newr):
		ii=np.where(r==r0)
		ii=ii[0][0]
		newmatrix[i,:]=matrix[ii,:]
	longdic={}
	for ri in np.arange(0,152,7.5):
		shortdicc={}
		for i,r0 in enumerate(newr):
			if r0-ri <7.5 and r0-ri>0:
				shortdicc[r0]=newmatrix[i,:]
			#else:
			#	print(r0)
		rlen=len(shortdicc.keys())
		if rlen==0:
			continue
		#print(rlen)
		A=np.zeros((rlen,newmatrix.shape[1]))
		for i,key in enumerate(shortdicc.keys()):
			A[i,:]=shortdicc[key]
		#	plt.plot(A[i,:],H,label=key)
		#plt.plot(np.nanmean(A,axis=0),H,label='mean')
		#plt.legend()
		#plt.show()
		longdic[ri]=np.nanmean(A,axis=0)
	AA=np.zeros((len(longdic.keys()),newmatrix.shape[1]))
#	print(AA.shape)
	rr=[]
	for j,key in enumerate(longdic.keys()):
		AA[j,:]=longdic[key]
		rr.append(key)
	#plt.contourf(AA.T,cmap='rainbow')
	#plt.show()
	return rr,AA
