"""

Gradient wind balance

"""


import os
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
from toolbox import dp_dr,getgradwind
from scipy.interpolate import griddata

storm=sys.argv[1]
sdate=sys.argv[2]
filename=sys.argv[3]
matrix=np.genfromtxt(filename)
r=matrix[:,0]
H=matrix[:,5]
u_cartesian=matrix[:,2]
v_cartesian=matrix[:,3]
v_tang=matrix[:,10]
x=matrix[:,12]
y=matrix[:,13]
fcor=matrix[:,9]
pressure=matrix[:,8]
thetas=matrix[:,1]
show=sys.argv[4]

if show == 'yes':
	show=True
else:
	show=False



ri=np.arange(0,120,5)
height_vec=Hi=np.arange(70,2500,50)
gradwind3d=np.zeros((len(ri)-1,len(height_vec)))
fcoriols2d=np.zeros((len(ri),len(height_vec)))
print('Getting radial gradient of pressure and gradent wind')
for jindex,hh in enumerate(height_vec):
  # Selection
    indices=np.where((H>=hh-25)& (H<hh+25))[0]
  # Quadratic interpolation
    interp_pres=np.polyfit(r[indices],pressure[indices], deg=3)

    a,b,c,d=interp_pres
    pres=a*ri**3+b*ri**2+c*ri+d
  # Linear interpolation of coriolis.

    finter = np.polyfit(r[indices], fcor[indices],1)
    fis=finter[0]*(ri)+finter[1]

  #Radial gradient of pressure
    dpdr,newradius=dp_dr(pres,ri)

  # Unit change of f
    fis=fis*10**(-2)

    fcoriols2d[:,jindex]=fis
  #Gradient wind function
    gradwind=getgradwind(dpdr,newradius,fis)
    gradwind3d[:,jindex]=gradwind

newradius=newradius/1000
plt.figure(figsize=(13,12))
ax=plt.subplot(311)
cs=ax.contourf(newradius,height_vec,gradwind3d.T,levels=np.arange(-5,np.nanmax(v_tang)-10,5),cmap='rainbow')
plt.colorbar(cs)
mean_azi=scipy.interpolate.griddata((r,H),v_tang, (newradius[None,:], Hi[:,None]),method='linear')
ax.set_title('Gradient ',fontsize=15)
ax.set_ylabel('Height [m]',fontsize=15)
ax=plt.subplot(312)
cs=ax.contourf(newradius,Hi,mean_azi,levels=np.arange(0,np.nanmax(v_tang)+2,5),cmap='jet')
plt.colorbar(cs)
ax.set_title('Azimuthal ',fontsize=15)
ax.set_ylabel('Height [m]',fontsize=15)
anomaly=mean_azi-gradwind3d.T
ax=plt.subplot(313)
cs=ax.contourf(newradius,Hi,anomaly,levels=np.arange(-np.nanmax(anomaly),np.nanmax(anomaly),4),cmap='bwr')
plt.colorbar(cs)
ax.set_title('Anomaly',fontsize=15)
ax.set_xlabel('Radius [km]',fontsize=15)
ax.set_ylabel('Height [m]',fontsize=16)
plt.suptitle('Gradiend wind balance '+storm+' on '+sdate,fontsize=18)
plt.savefig('figs/'+storm+'/axisym/gradwind'+sdate+'.png')
if show:
	plt.show()
plt.close()


momentum=np.zeros((len(ri),len(height_vec)))
vt_interp=scipy.interpolate.griddata((r,H),v_tang, (ri[None,:], Hi[:,None]),method='linear')

#Momentum estimation
for index,r0 in enumerate(ri):
    vr=vt_interp[:,index]*r0*1000
    fr2=fcoriols2d[index,:]*(r0**2)
    momentum[index,:]=vr+(0.5)*fr2

plt.contourf(ri,Hi,momentum.T,levels=10**(6)*np.arange(0,7,1),cmap='Spectral')
plt.colorbar(label=r'$m^2 s^{-1}$')
plt.ylabel('Height [m]',fontsize=15)
plt.xlabel('Radius [km]',fontsize=15)
plt.title('Momentum on '+storm+' on '+sdate,fontsize=16 )
plt.savefig('figs/'+storm+'/axisym/momentum'+sdate+'.png')
if show:
	plt.show()


quit()
