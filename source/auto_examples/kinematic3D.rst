

.. _sphx_glr_auto_examples_kinematic3D.py:


Plot the 3D kinematic structure of a TC
=====================================================

One of the most important kinematic metrics are vorticity and divergence. Both
directly related to the first spatial derivative of the wind components.






.. code-block:: python


    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from toolbox import convergence_vorticity
    from scipy.interpolate import griddata

    matrix=np.genfromtxt('tempjulia.txt')
    radius=matrix[:,0]
    Height=H=matrix[:,5]
    u_cartesian=matrix[:,2]
    v_cartesian=matrix[:,3]
    u_radial=matrix[:,11]
    v_tang=matrix[:,10]
    x=matrix[:,12]
    y=matrix[:,13]
    temperature=matrix[:,6]
    pressure=matrix[:,8]
    thetas=matrix[:,1]

    storm=sys.argv[1]
    sdate=sys.argv[2]
    print(sdate)
    figdir='/home/jlgf/Documents/MRes/Project/figs/'+storm+'/'
    xi=np.arange(-140,140,5)
    yi=np.arange(-140,140,5)
    height_vec=Hi=np.arange(50,2800,50)
    vorticity=np.zeros((len(xi),len(yi),len(Hi)))
    divergence=np.zeros((len(xi),len(yi),len(Hi)))


    for ij,H0 in enumerate(height_vec):
        try:
            indices=np.where((H>=H0-50)& (H<H0))[0]
            u_xh=griddata((x[indices],y[indices]),u_cartesian[indices], (xi[None,:], yi[:,None]),method='linear')
            v_xh=griddata((x[indices],y[indices]),v_cartesian[indices], (xi[None,:], yi[:,None]),method='linear')
        except:
            vorticity[:,:,ij]=divergence[:,:,ij]=np.nan
            continue
        nabla_dot_u,nabla_cross_u=convergence_vorticity(u_xh,v_xh,xi,yi,dim_order='xy')
        vorticity[:,:,ij]=nabla_cross_u*100
        divergence[:,:,ij]=nabla_dot_u*100


    plt.figure(figsize=(14,10))
    ax=plt.subplot(221)
    CS=ax.contourf(xi,Hi,vorticity[:,np.where(yi==0)[0][0],:].T,cmap='bwr',levels=np.linspace(np.nanmin(vorticity),np.nanmax(vorticity),10))
    #plt.colorbar(CS)
    plt.xlim([-45,45])
    plt.ylim([0,2800])
    plt.title('Vertical vorticity W-E Height Cross section')
    ax=plt.subplot(222)
    CS=ax.contourf(yi,Hi,vorticity[np.where(xi==0)[0][0],:,:].T,cmap='bwr',levels=np.linspace(np.nanmin(vorticity),np.nanmax(vorticity),10))
    plt.colorbar(CS,label=r'$10^{-2}$s$^{-1}$')
    plt.xlim([-45,45])
    plt.ylim([0,2800])
    plt.title('Vertical vorticity N-S Height Cross section')
    ax=plt.subplot(223)
    CS=ax.contourf(yi,Hi,divergence[:,np.where(yi==0)[0][0],:].T,cmap='bwr',levels=np.linspace(np.nanmin(divergence),np.nanmax(divergence),10))
    #plt.colorbar(CS)
    plt.xlim([-45,45])
    plt.ylim([0,2800])
    plt.title('Horizontal divergence W-E Height Cross section')
    ax=plt.subplot(224)
    CS=ax.contourf(yi,Hi,divergence[np.where(xi==0)[0][0],:,:].T,cmap='bwr',levels=np.linspace(np.nanmin(divergence),np.nanmax(divergence),10))
    plt.colorbar(CS,label=r'$10^{-2}$s$^{-1}$')
    plt.xlim([-45,45])
    plt.ylim([0,2800])
    plt.title('Horizontal divergence N-S Height Cross section')
    plt.suptitle('Kinematic plots for '+storm+' on '+sdate,fontsize=15)
    plt.savefig(figdir+'axisym/kinematic'+sdate+'.png')
    #plt.close()
    plt.show()
    # Radial average
    ri=np.arange(0,120,10)
    #u_xh=scipy.interpolate.griddata((radius,Height),u_cartesian[indices], (ri[None,:], Hi[:,None]),method='linear')
    #v_xh=scipy.interpolate.griddata((radius,Height),v_cartesian[indices], (ri[None,:], Hi[:,None]),method='linear')

    # Plan views at 150 m, 400 m, 800 m, 2000 m,
    Heights=[150,400,800,2000]
    fig=plt.figure(figsize=(18,8))
    for index,height in enumerate(Heights):
        ax=plt.subplot(241+index)
        try:
            CS=ax.contourf(xi,yi,vorticity[:,:,np.where(Hi==height)[0][0]],cmap='bwr',levels=np.linspace(np.nanmin(vorticity[:,:,np.where(Hi==height)[0][0]]),np.nanmax(vorticity),12))
        except:
            continue
        #plt.colorbar(CS)
        plt.colorbar(CS,label=r'$10^{-2}$s$^{-1}$')
        plt.xlim([-45,45])
        plt.ylim([-45,45])
        plt.title('Vertical vorticity at '+str(height)+' m')
        ax=plt.subplot(245+index)
        CS=ax.contourf(xi,yi,divergence[:,:,np.where(Hi==height)[0][0]],cmap='bwr',levels=np.linspace(np.nanmin(divergence[:,:,np.where(Hi==height)[0][0]]),np.nanmax(divergence),12))
        plt.colorbar(CS,label=r's$^{-1}$')
        plt.xlim([-45,45])
        plt.ylim([-45,45])
        plt.title('Horizontal divergence at '+str(height)+' m')
    plt.suptitle('Kinematic plots for '+storm+' on '+sdate,fontsize=15)#
    plt.savefig(figdir+'planviews/plankinematic'+sdate+'.png')
    plt.close()

**Total running time of the script:** ( 0 minutes  0.000 seconds)



.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: kinematic3D.py <kinematic3D.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: kinematic3D.ipynb <kinematic3D.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.readthedocs.io>`_
