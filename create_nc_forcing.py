! create nc file for large scale SCM forcing

import numpy as np
from netCDF4 import Dataset
import datetime
from datetime import date, timedelta

###### Original forcing file
file="/Users/uma2103/Downloads/mao180varanaecmwfanaradarC1.c1.20140201.000000.cdf"
nc=Dataset(file);


timestart=256; #corrosponds to March 19, hour=0
timeend=274;   #corrosponds to March 21, hour=6 [actully 271 is March 20, hour=21]


########## TIMES #############
TIME=[];
DT=datetime.datetime(2014, 3, 19, 0, 0, 0)
tim=DT.strftime("%Y-%m-%d_%H:%M:%S")
TIM=[tim[0],tim[1],tim[2],tim[3],tim[4],tim[5],tim[6],tim[7],tim[8],tim[9],tim[10],tim[11],tim[12],tim[13],tim[14],tim[15],tim[16],tim[17],tim[18]]

TIME.append(TIM)

for i in range(0,17):
    DT=DT+timedelta(hours=3)
    tim=DT.strftime("%Y-%m-%d_%H:%M:%S")
    TIM=[tim[0],tim[1],tim[2],tim[3],tim[4],tim[5],tim[6],tim[7],tim[8],tim[9],tim[10],tim[11],tim[12],tim[13],tim[14],tim[15],tim[16],tim[17],tim[18]]
    TIME.append(TIM)

TIME_F=np.array(TIME)

###############################
# calculation of rho
rgas=287.
grav=9.81

rho=np.zeros([18,45])
h=np.zeros([18,45]) #height


p=np.array(nc.variables['lev'])
T=np.array(nc.variables['T'])
TS=np.array(nc.variables['T_skin'])

for i,k in enumerate(range(timestart,timeend)):
    rho[i,:]=p*100/T[k,:]/rgas
    h[i,:]=66 + (TS[i]+273.15)/-0.0065 *( (p/p[0])**(-(8.31431*-0.0065/9.806/0.0289644)) -1)

    #z[i,:] = -1/rho[i,:]/grav*-25*100 # -25 is dp

#create nc file
rootgrp = Dataset("/Users/uma2103/GoAmazon_March19_March20_FORCE.nc", "w", format="NETCDF4")

level = rootgrp.createDimension("force_layers", 45)
time = rootgrp.createDimension("Time", None)
DateStrLen= rootgrp.createDimension("DateStrLen", None) # 18 time step

times = rootgrp.createVariable("time","f8",("Time",))
levels = rootgrp.createVariable("level","i4",("force_layers",))


## Times
Times = rootgrp.createVariable("Times","S1",("Time","DateStrLen"))

# fill in 
Times[:,:]= TIME_F[:,:];


# Z
Z_FORCE = rootgrp.createVariable("Z_FORCE","f4",("Time","force_layers"))
Z_FORCE.FieldType=" 104.0"
Z_FORCE.MemoryOrder=" Z"
Z_FORCE.description=" height of forcing time series"
Z_FORCE.units=" m "
# fill in 
Z_FORCE[:,:]= h[:,:]


W_SUBS = rootgrp.createVariable("W_SUBS","f4",("Time","force_layers"))
W_SUBS.FieldType=" 104.0"
W_SUBS.MemoryOrder=" Z"
W_SUBS.description=" large-scale vertical velocity"
W_SUBS.units=" m s-1"
# fill in 
WW_SUBS=np.zeros([18,45])
WW_SUBS[:,1:45]= np.array(nc.variables['omega'])[timestart:timeend,0:44]*(100/3600.)

W_SUBS[:,:]= -WW_SUBS/rho/grav #this is now m/s


#theta large scale (or dry static energy?)
TH_LARGESCALE = rootgrp.createVariable("TH_LARGESCALE","f4",("Time","force_layers"))
TH_LARGESCALE.FieldType=" 104.0"
TH_LARGESCALE.MemoryOrder=" Z"
TH_LARGESCALE.description=" theta large-scale"
TH_LARGESCALE.units=" K"
# fill in 
TH_LARGESCALE[:,:]= np.array(nc.variables['s'])[timestart:timeend,:]



TH_LARGESCALE_TEND = rootgrp.createVariable("TH_LARGESCALE_TEND","f4",("Time","force_layers"))
TH_LARGESCALE_TEND.FieldType=" 104.0"
TH_LARGESCALE_TEND.MemoryOrder=" Z"
TH_LARGESCALE_TEND.description=" horizontal and vertical temperature tendency"
TH_LARGESCALE_TEND.units=" K s-1"
# fill in 
TH_LARGESCALE_TEND[:,:]= -(np.array(nc.variables['s_adv_h'][timestart:timeend,:])+np.array(nc.variables['s_adv_v'][timestart:timeend,:]))/(3600.)


QV_LARGESCALE = rootgrp.createVariable("QV_LARGESCALE","f4",("Time","force_layers"))
QV_LARGESCALE.FieldType=" 104.0"
QV_LARGESCALE.MemoryOrder=" Z"
QV_LARGESCALE.description=" Qv"
QV_LARGESCALE.units=" kg kg-1"
# fill in 
QV_LARGESCALE[:,:]= (np.array(nc.variables['q'])[timestart:timeend,:])/1000.


QV_LARGESCALE_TEND = rootgrp.createVariable("QV_LARGESCALE_TEND","f4",("Time","force_layers"))
QV_LARGESCALE_TEND.FieldType=" 104.0"
QV_LARGESCALE_TEND.MemoryOrder=" Z"
QV_LARGESCALE_TEND.description=" large scale qv tendency"
QV_LARGESCALE_TEND.units=" kg kg-1 s-1"
# fill in 
QV_LARGESCALE_TEND[:,:]= -(np.array(nc.variables['q_adv_h'][timestart:timeend,:])+np.array(nc.variables['q_adv_v'][timestart:timeend,:]))/(3600.)



U_LARGESCALE = rootgrp.createVariable("U_LARGESCALE","f4",("Time","force_layers"))
U_LARGESCALE.FieldType=" 104.0"
U_LARGESCALE.MemoryOrder=" Z"
U_LARGESCALE.description=" U"
U_LARGESCALE.units=" m s-1"
# fill in 
U_LARGESCALE[:,:]= np.array(nc.variables['u'])[timestart:timeend,:]


V_LARGESCALE = rootgrp.createVariable("V_LARGESCALE","f4",("Time","force_layers"))
V_LARGESCALE.FieldType=" 104.0"
V_LARGESCALE.MemoryOrder=" Z"
V_LARGESCALE.description=" V"
V_LARGESCALE.units=" m s-1"
# fill in 
V_LARGESCALE[:,:]= np.array(nc.variables['v'])[timestart:timeend,:]



TAU_LARGESCALE = rootgrp.createVariable("TAU_LARGESCALE","f4",("Time","force_layers"))
TAU_LARGESCALE.FieldType=" 104.0"
TAU_LARGESCALE.MemoryOrder=" Z"
TAU_LARGESCALE.description=" largescale timescale"
TAU_LARGESCALE.units=" s"
# fill in 
TAU_LARGESCALE[:,:]= 3600*np.ones([18,45])



TAU_LARGESCALE_TEND = rootgrp.createVariable("TAU_LARGESCALE_TEND","f4",("Time","force_layers"))
TAU_LARGESCALE_TEND.FieldType=" 104.0"
TAU_LARGESCALE_TEND.MemoryOrder=" Z"
TAU_LARGESCALE_TEND.description=" tendency largescale timescale"
TAU_LARGESCALE_TEND.units=" s"
# fill in 
TAU_LARGESCALE_TEND[:,:]= np.zeros([18,45])


HFX_FORCE = rootgrp.createVariable("HFX_FORCE","f4",("Time"))
HFX_FORCE.FieldType=" 104.0"
HFX_FORCE.MemoryOrder=" Z"
HFX_FORCE.description=" SCM ideal surface sensible heat flux"
HFX_FORCE.units=" W m-2"
# fill in 
HFX_FORCE[:]= np.array(nc.variables['SH'])[timestart:timeend]


LH_FORCE = rootgrp.createVariable("LH_FORCE","f4",("Time"))
LH_FORCE.FieldType=" 104.0"
LH_FORCE.MemoryOrder=" Z"
LH_FORCE.description=" SCM ideal surface latent heat flux"
LH_FORCE.units=" W m-2"
# fill in 
LH_FORCE[:]= np.array(nc.variables['LH'])[timestart:timeend]



TSK_FORCE = rootgrp.createVariable("TSK_FORCE","f4",("Time"))
TSK_FORCE.FieldType=" 104.0"
TSK_FORCE.MemoryOrder=" T"
TSK_FORCE.description=" SCM ideal surface skin temperature"
TSK_FORCE.units=" K"
# fill in 
TSK_FORCE[:]= np.array(nc.variables['T_skin'])[timestart:timeend]+273.15


rootgrp.close()
######### set time: March 19 - March 20 2014 [256:274]


