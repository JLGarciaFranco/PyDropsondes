import os
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
from toolbox import dp_dr,getgradwind
from scipy.interpolate import griddata,RegularGridInterpolator
import metpy.calc as mpcalc
from metpy.units import units

storm=sys.argv[1]
sdate=sys.argv[2]
figdir='/home/jlgf/Documents/MRes/Project/figs/'+storm+'/'
matrix=np.genfromtxt('tempjulia.txt')
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


ri=np.arange(0,120,5)
height_vec=Hi=np.arange(50,2500,50)
gradwind3d=np.zeros((len(ri)-1,len(height_vec)))
fcoriols2d=np.zeros((20,len(height_vec)))

print('Getting radial gradient of pressure and gradent wind')
for jindex,hh in enumerate(height_vec):
    indices=np.where((H>=hh-25)& (H<hh+25))[0]
    interp_pres=np.polyfit(r[indices],pressure[indices], deg=2)
    a,b,c=interp_pres
    pres=a*ri**2+b*ri+c


    finter = np.polyfit(r[indices], fcor[indices],1)
    fis=finter[0]*(ri)+finter[1]
    #Radial gradient of pressure
    dpdr,newradius=dp_dr(pres,ri)

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
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') #
mean_azi=scipy.interpolate.griddata((r,H),v_tang, (newradius[None,:], Hi[:,None]),method='linear')
secondazi=scipy.interpolate.griddata((r,H),v_tang, (ri[None,:], Hi[:,None]),method='linear')
ax.set_ylabel('Height [m]')
ax.set_title('Gradient ')
ax=plt.subplot(312)
cs=ax.contourf(newradius,Hi,mean_azi,levels=np.arange(0,np.nanmax(v_tang)+2,5),cmap='jet')
plt.colorbar(cs)
ax.set_title('Azimuthal ')
ax.set_ylabel('Height [m]')
anomaly=mean_azi-gradwind3d.T
ax=plt.subplot(313)
cs=ax.contourf(newradius,Hi,anomaly,levels=np.arange(-np.nanmax(anomaly),np.nanmax(anomaly),4),cmap='bwr')
plt.colorbar(cs)
ax.set_title('Anomaly')
ax.set_xlabel('Radius [km]')
ax.set_ylabel('Height [m]')
plt.suptitle('Gradiend wind balance '+storm+' on '+sdate,fontsize=16)
#plt.savefig(figdir+'axisym/gradwind'+sdate+'.png')
#plt.show()
plt.close()
momentum=np.zeros((len(ri),len(height_vec)))
#Momentum estimation
print(len(ri))
for index,r0 in enumerate(newradius):
    vr=secondazi[:,index]*r0*1000
    #x=yy
    fr2=fcoriols2d[index,:]*(r0**2)
    momentum[index,:]=vr+(0.5)*fr2
plt.contourf(ri,Hi,momentum.T,levels=10**(6)*np.arange(0,7,1),cmap='Spectral')
plt.colorbar()
plt.ylabel('Height [m]',fontsize=15)
plt.xlabel('Radius [km]',fontsize=15)
plt.title('Momentum on '+storm+' on '+sdate,fontsize=16 )
plt.savefig(figdir+'axisym/momentum'+sdate+'.png')
plt.show()


quit()




plt.scatter(r[indices],v_tang[indices],color='red',label=r'$V_t$')
plt.plot(newradius/1000.,gradwind,'--k',label=r'$V_g$')
plt.grid(alpha=0.5)
plt.xlabel('Radius [km]')
plt.ylabel(r'Wind speed m s$^{-1}$')
plt.suptitle('Gradiend wind balance Earl')
plt.legend()
plt.show()
#        continue
plt.scatter(r[indices],pressure[indices],color='red',label='Pressure Obs.')
plt.plot(ri,pres,'--k',label='Squared fit')
plt.grid(alpha=0.5)
plt.xlabel('Radius [km]')
plt.ylabel(r'Pressure $[mb]$')
plt.title('Pressure Earl 28 August')
plt.legend()
plt.show()
