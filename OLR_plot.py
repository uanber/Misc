#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import colorConverter
from mpl_toolkits.basemap import Basemap, cm, interp as BasemapInterp
import numpy as np
from netCDF4 import Dataset
from netcdftime import utime
from datetime   import datetime
import os, errno

lonrange = [250., 5]
latrange = [-2.5,   38.]

runname='C384n3-GHS_fore_nh_0313_050701_oneway'
plotname='slp'
dir = '/archive/Usama.Anber/Hiram/C384n3-GHS_fore_nh_0313_050701_oneway/history/20050802/'
coarsedata = dir + '20050802.atmos_8xdaily.nc'
#coarsegrid = dir + '20050901.grid_spec.tile6.nc'
nesteddata = dir + '20050802.atmos_8xdaily.nest02.nc'
nestedgrid = dir + '20050802.grid_spec.nest02.nc'

fcd = Dataset(coarsedata,'r')
#fcg = Dataset(coarsegrid,'r')

# grid_lon  = fcg.variables['grid_lon'][:]
# grid_lat  = fcg.variables['grid_lat'][:]
gxt  = fcd.variables['grid_xt'][:]
gyt  = fcd.variables['grid_yt'][:]
grid_lon, grid_lat = np.meshgrid(gxt,gyt)
olr   = fcd.variables['slp'][:]

times = fcd.variables['time'][:]
cdftime = utime(getattr(fcd.variables['time'],'units'))

fnd = Dataset(nesteddata,'r')
fng = Dataset(nestedgrid,'r')

grid_lon_n  = fng.variables['grid_lon'][:,:]
grid_lat_n  = fng.variables['grid_lat'][:,:]
olr_n   = fnd.variables['slp'][:]

#Create alpha colormap
# http://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
color1 = colorConverter.to_rgba('white')
color2 = colorConverter.to_rgba('black')
cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',[color1,color2],256)
cmap._init()
alphas = np.flipud(np.linspace(0,1,cmap.N+3))
cmap._lut[:,-1] = alphas

for t,time in enumerate(times):

    datestr = cdftime.num2date(time)
    print datestr

    figM = plt.figure()
    mM = Basemap(llcrnrlon=lonrange[0],   llcrnrlat=latrange[0],
                 urcrnrlon=lonrange[-1],  urcrnrlat=latrange[-1],
                 projection='gnom', lat_0=22.5,lon_0=-60,
                 resolution='i')
    mM.bluemarble()
    mM.drawcoastlines()
    mM.drawcountries()

    if (t == 0):
       x,y = mM(grid_lon[:,:],grid_lat[:,:])

    IR  = olr[t,:,:]
    IR = np.ma.masked_where(IR >= 200, IR)   
    im = mM.pcolormesh(x,y,IR,cmap=plt.cm.Greys,vmin=80,vmax=220,edgecolors='none',shading='flat')

    #mM.drawcoastlines()
    #mM.drawcountries()
    #parallels = np.arange(0.,80,10.)
    #mM.drawparallels(parallels,labels=[0,0,0,0])
    #meridians = np.arange(10.,360.,10.)
    #mM.drawmeridians(meridians,labels=[0,0,0,0])
    #plt.colorbar()
    plt.title(datestr,loc='left',fontsize=10,)

    plt.savefig('coarsegrid-' + runname + '-' + plotname  + '-%03d.png' % (t+1),dpi=150, bbox_inches='tight', transparent=True)

    if (t == 0): 
        x_n,y_n = mM(grid_lon_n[:,:],grid_lat_n[:,:])

    IR_n = olr_n[t,:,:]
    IR_n = np.ma.masked_where(IR_n >= 200., IR_n)
    im = mM.pcolormesh(x_n,y_n,IR_n,cmap=plt.cm.Greys,vmin=80,vmax=220,shading='flat')
    mM.plot(x_n[0,:],y_n[0,:],'pink',linewidth=3)
    mM.plot(x_n[-1,:],y_n[-1,:],'pink',linewidth=3)
    mM.plot(x_n[:,0],y_n[:,0],'pink',linewidth=3)
    mM.plot(x_n[:,-1],y_n[:,-1],'pink',linewidth=3)

    #mM.drawparallels(parallels,labels=[0,0,0,0])
    #mM.drawmeridians(meridians,labels=[0,0,0,0])

    plt.savefig('nestedgrid-' + runname + '-' + plotname  + '-%03d.png' % (t+1),dpi=150, bbox_inches='tight')
    plt.show()
    plt.close('all')
