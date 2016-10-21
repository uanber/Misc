# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:20:39 2016

@author: uma2103
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:38:29 2016

@author: uma2103
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import date, timedelta
import datetime
BASIN='EP'#'NA'
STR='HU'#'HU'#'TS'

RANGE_YEARS=range(2005,2015)
F=[];
for YEAR in RANGE_YEARS:
    file="/Users/uma2103/IBTr/IBTr/"+str(YEAR)+"_IBTr_wmo_"+BASIN+STR+".txt"

    df=pd.read_csv(file, header=None, skiprows=[0], sep=r"\s+")
    df.columns=['date','sfc wind','lon','lat','Ps','NAN','dont_know', 'name']

    STORM_NAME=[]
    for i in range(0,len(df['name'])):
        if (i != len(df['name'])-1 and df['name'][i] !=df['name'][i+1]):
            STORM_NAME.append(df['name'][i])
        elif (i==len(df['name'])-1):  
            STORM_NAME.append(df['name'][i])
        

        FORE_D1=datetime.date(int(YEAR), 7, 1) #1st date simulations
        FORE_D2=FORE_D1+timedelta(days=14)
        TARGET_DATE=datetime.date(int(YEAR), 11, 4) #last date simulations
        
        S=[] # array for mean number of storms in each lead
        while FORE_D2 <= TARGET_DATE:
          

            MASK=[df['name']==i for i in STORM_NAME]

            while (FORE_D2 <= TARGET_DATE):
                C=0.
        
                for kk in range(len(MASK)): #loop over the entire storms

                    STORM_DATE= np.array(df[MASK[kk]]['date'])

                    FIRST=str(STORM_DATE[0]) # storm first date
                    LAST=str(STORM_DATE[-1]) # storm last date
                    # now convert date from string to datetime 
                    STORM_FIRST_DATE=datetime.date(year=int(FIRST[0:4]), month=int(FIRST[4:6]), day=int(FIRST[6:8]))
                    STORM_LAST_DATE=datetime.date(year=int(LAST[0:4]), month=int(LAST[4:6]), day=int(LAST[6:8]))

                    #print "first ", STORM_FIRST_DATE
                    #print "last ", STORM_LAST_DATE
        
                    if (FORE_D1<STORM_FIRST_DATE<FORE_D2 or FORE_D1<STORM_LAST_DATE<FORE_D2):
                        C=C+1
                #print "Number of obsereved storms between: ", FORE_D1.strftime('%m-%d'), " and ", FORE_D2.strftime('%m-%d'), " is: ", C
                FORE_D1=FORE_D1+timedelta(days=5)
                FORE_D2=FORE_D1+timedelta(days=14)
                S.append(C)
    F.append(S)

# plotting 10 years of sub-seasonal runs.
M=np.zeros(240)

for i in range(0,10):
    M[24*i:23+24*i]=F[i]
    
plt.plot(M)    