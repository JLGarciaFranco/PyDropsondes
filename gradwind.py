import numpy
import matplotlib.pyplot as plt
import os
import sys
def gradientwind(pressure,v_azimuth,rdistance):
    plt.figure(figsize=(12,12))
    ax=plt.subplot(211)
    ax.scatter(rdistance,pressure,marker='v',color='black')
    ax.set_xlabel('Distance [km]',fontsize=12)
    ax.set_ylabel('Pressure [hPa]',fontsize=12)
    ax.grid()
    ax=plt.subplot(212)
    ax.scatter(rdistance,v_azimuth,marker='v',color='black')
    ax.set_xlabel('Distance [km]',fontsize=12)
    ax.set_ylabel(r'Azimuthal wind speed [$ms^{-1}$]',fontsize=12)
    ax.grid()
    plt.suptitle('Isabel 09/12/03 Pressure analysis',fontsize=14)
    plt.savefig('/home/jlgf/Documents/MRes/Project/figs/Isabel/gradwind.png')
    plt.show()
