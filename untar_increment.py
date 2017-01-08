#!/usr/bin/env python

import os
import tarfile
from datetime import datetime, timedelta
from datetime import date
import datetime

year='2013'

dd= datetime.date(2013, 7, 1)
five_day_increment = datetime.timedelta(days=5)
half_month_increment = datetime.timedelta(days=15)

while (dd< datetime.date(2013, 9, 4)):
#while (dd< datetime.date(2013, 10, 9)):

#for i in range(0,len(dd)):
    kk = dd.strftime('%m%d')
    print kk
    #m = dd + half_month_increment # the 15 or 5 day increment 
    #mm = m.strftime('%m%d') #convert date to string
    #print 'folder_',kk,'_file_',mm
    
    #folder_path=("/archive/Usama.Anber/Hiram/C384n3-GHS_subfore_nh_0313_"+year[2:4]+kk+"-C1/history/"+year+mm)
    folder_path = ("/archive/Usama.Anber/Hiram/C384n3-GHS_subfore_nh_0313_"+year[2:4]+kk+"-C3/history/"+year+kk)
    print folder_path
    os.makedirs(folder_path)   # make directory
    #tar_file=("/archive/Usama.Anber/Hiram/C384n3-GHS_subfore_nh_0313_"+year[2:4]+kk+"-C1/history/"+year+mm+".nc.tar")     # name of the tar file
    tar_file    = ("/archive/Usama.Anber/Hiram/C384n3-GHS_subfore_nh_0313_"+year[2:4]+kk+"-C3/history/"+year+kk+".nc.tar")
    tar=tarfile.open(tar_file)
    tar.extractall(folder_path)  # extact the files into the file_path folder
    tar.close()

    dd = dd + five_day_increment
