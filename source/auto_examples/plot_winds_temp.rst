

.. _sphx_glr_auto_examples_plot_winds_temp.py:


Plotting the wind and temperature structure of a TC
=====================================================

This example gives the typical plot of pressure and wind speed from the best track dataset.

Specifically, the `best track <https://www.nhc.noaa.gov/data/#hurdat>`_ dataset, provides the maximum sustained winds at 10 m altitude (hereafter :math:`U_{10}`)
and the minimum surface pressure (hereafter :math:`P_{min}`).





.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/images/sphx_glr_plot_winds_temp_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/images/sphx_glr_plot_winds_temp_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/images/sphx_glr_plot_winds_temp_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out::

    21-09-2005




|


.. code-block:: python

    import os
    import numpy as np
    import scipy
    import sys
    import matplotlib.pyplot as plt
    from toolbox import potential_temperature
    from scipy.interpolate import griddata
    matrix=np.genfromtxt('tempjulia.txt')
    r=matrix[:,0]
    H=matrix[:,5]
    u_radial=matrix[:,11]
    v_tang=matrix[:,10]
    x=matrix[:,12]
    y=matrix[:,13]
    temperature=matrix[:,6]
    thetae=matrix[:,7]
    binning=False
    thetas=matrix[:,1]
    #storm=sys.argv[1]
    #sdate=sys.argv[2]
    storm='Rita'
    sdate='21-09-2005'
    print(sdate)
    figdir='/home/jlgf/Documents/MRes/Project/figs/'+storm+'/'

    ri=np.arange(0,100,2.5)
    Hi=np.arange(50,3200,100)
    mean_azi=scipy.interpolate.griddata((r,H),v_tang, (ri[None,:], Hi[:,None]),method='linear')
    mean_radial=scipy.interpolate.griddata((r,H),u_radial, (ri[None,:], Hi[:,None]),method='linear')
    mean_temp=scipy.interpolate.griddata((r,H),temperature, (ri[None,:], Hi[:,None]),method='linear')
    mean_vert=scipy.interpolate.griddata((r,H),thetae, (ri[None,:], Hi[:,None]),method='linear')
    plt.figure(figsize=(14,10))
    plottingdictionary={"Radial wind":mean_radial,"Azimuthal wind":mean_azi,"Potential temperature":mean_temp,r"$\theta_e$":mean_vert}
    colormaps=['seismic','rainbow','coolwarm','gist_rainbow']
    labels=[r'm s$^{-1}$',r'm s$^{-1}$','K',r'K']
    spacing=[2,4,1.5,2]
    counter=0
    for variable in plottingdictionary.keys():
    	field=plottingdictionary[variable]
    	ax=plt.subplot(221+counter)
    	CS=plt.contourf(ri,Hi,field,cmap=colormaps[counter],levels=np.arange(int(np.nanmin(field))-(spacing[counter]/2),int(np.nanmax(field))+(spacing[counter]),spacing[counter]))
    	plt.xlabel('Radius [km] ',fontsize=14)
    	plt.ylabel('Height [m] ',fontsize=14)
    #	ax.scatter(r,H,color='black',s=1)
    	plt.title(variable,fontsize=16)
    	plt.xlim([0,100])
    	plt.ylim([0,3000])
    	plt.colorbar(CS,label=labels[counter])
    	counter+=1

    plt.suptitle(' Cross sections of '+storm+' on '+str(sdate),fontsize=19)
    plt.savefig('figs/crossect_nolines'+str(sdate)+'.png')

    #	print(indices[0])
    heights=[100,400,800,2000]
    figwinds=plt.figure(figsize=(18,7))
    figtemps=plt.figure(figsize=(18,7))
    xi=np.arange(-95,95,2.5)
    yi=np.arange(-95,95,2.5)
    for counter,hh in enumerate(heights):
    	indices=np.where((H>hh-40)& (H<hh+40))
    	shortx=x[indices]
    	shorty=y[indices]
    	spacing=[2,4,0.75,1.5]
    	plottingdictionary={"Radial wind":u_radial,"Azimuthal wind":v_tang,"Potential temperature":temperature,r"$\theta_e$":thetae}
    	for incounter,key in enumerate(plottingdictionary.keys()):
    		if key=="Radial wind":
    			ax=figwinds.add_subplot(241+counter)
    		elif key == "Azimuthal wind":
    			ax=figwinds.add_subplot(245+counter)
    		elif key=="Potential temperature":
    			ax=figtemps.add_subplot(241+counter)
    		else:
    			ax=figtemps.add_subplot(245+counter)
    		sliced_var=plottingdictionary[key][indices]
    		interp_var=scipy.interpolate.griddata((shortx,shorty),sliced_var, (xi[None,:], yi[:,None]),method='linear')
    		if key=="Radial wind" and hh <500:
    			levelss=np.arange(int(np.nanmin(interp_var))-(spacing[incounter]/2),-int(np.nanmin(interp_var))+spacing[incounter],spacing[incounter])
    		else:
    			levelss=np.arange(int(np.nanmin(interp_var)-(spacing[incounter]/2)),int(np.nanmax(interp_var)+(spacing[incounter])),spacing[incounter])
    		cs=ax.contourf(xi,yi,interp_var,cmap=colormaps[incounter],levels=levelss)
    		plt.colorbar(cs,ax=ax)
    		ax.set_title(key+' at '+str(hh)+' m',fontsize=16)
    figwinds.suptitle('Plan views of cylindrical winds of '+storm+' on '+str(sdate),fontsize=19)
    figtemps.suptitle('Plan views of temperature fields of '+storm+' on '+str(sdate),fontsize=19)
    figwinds.savefig('figs/winds_'+str(sdate)+'.png')
    figtemps.savefig('figs/temps_'+str(sdate)+'.png')

**Total running time of the script:** ( 0 minutes  2.843 seconds)



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: plot_winds_temp.py <plot_winds_temp.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: plot_winds_temp.ipynb <plot_winds_temp.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.readthedocs.io>`_
