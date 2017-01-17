import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from netCDF4 import Dataset

path='/archive/Usama.Anber/verona/ap-dp/HiRAM4-crm-2km-OML50///history/'
#HiRAM4-crm-2km-OML50/
#HiRAM4-crm-001km-OML
#date=('19800921','19801021','19801121', '19801221', '19810121', '19810221', '19810321', '19810421', '19810521', '19810621', '19810721', '19810821', '19810921', '19811021', '19811121', '19811221','19820121', '19820221', '19820321', '19820421', '19820521', '19820621', '19820721', '19820821', '19820921', '19821021', '19821121', '19821221', '19830121', '19830221')


date=('19800921','19801021','19801121', '19801221', '19810121', '19810221', '19810321', '19810421', '19810521', '19810621', '19810721', '19810821', '19810921', '19811021', '19811121', '19811221','19820121', '19820221')
LW=[];

for d in date:
    file=path+d+'/'+d+'.atmos_4xdaily.nc'
    print file
    nc=Dataset(file)
    lwdn_sfc=nc.variables['lwdn_sfc']
    lwdn_sfcM=np.mean(np.mean(lwdn_sfc,axis=2),axis=1)
    LW=LW+lwdn_sfcM.tolist()
    
    
# save LW to a txt file   
L=open('/home/Usama.Anber/LWUP_50.txt','w')
for ele in LW:
    L.write(str(ele)+'\n')

L.close()

plt.plot(LW); plt.show    

