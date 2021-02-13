import pandas as pd
import numpy as np
from scipy import stats
import math
import itertools
import time

input_data = 'neutral_200.csv'
input_mns_sds = 'input.xls'
#import data
og_data = pd.read_csv(input_data)
#import means and sds
input = pd.read_excel(input_mns_sds)
means_sds = {}
generated_avgs = {}
#gather means and sds
for month in range(12):
    month+=1
    #each month will have a weather dict with 8 keys for each weather variable (rain, min temp, max temp, srad)
    weather_dict = {'rain':[input.loc[month-1,'rain1']/input.loc[month-1,'rnum'],input.loc[month-1,'rain2']], 'tmin':[input.loc[month-1,'tmin'],input.loc[month-1,'tmin2']], 'tmax':[input.loc[month-1,'tmax'],input.loc[month-1,'tmax2']], 'srad':[input.loc[month-1,'srad'],input.loc[month-1,'srad2']], 'tnw':[input.loc[month-1,'tnmnw'],input.loc[month-1,'tnsdw']],'tnd':[input.loc[month-1,'tnmnd'],input.loc[month-1,'tnsdd']],'txw':[input.loc[month-1,'txmnw'],input.loc[month-1,'txsdw']],'txd':[input.loc[month-1,'txmnd'],input.loc[month-1,'txsdd']]}
    means_sds[month] = weather_dict
    generated_avgs[month] = {'rain':0,'tnw':0,'tnd':0,'txw':0,'txd':0}

#generate the parameters
b = 18.77942**2/means_sds[1]['rain'][0]
a = means_sds[1]['rain'][0]/b

#inputs a singular column of rainfall values, returns the probability of wet day given dry and prob of dry given wet
def gen_monthly_probs(rain_df):

    #creates an indicator column for if it rained that day or not
    ind = rain_df.apply(lambda x: [1 if y >0.01 else 0 for y in x]).reset_index().drop(columns = ['index'])

    #find the number of wet and dry days
    last_day = ind['pred'].iloc[-1]
    if last_day==1:
        num_wet = ind.sum(axis=0).values[0]-1
        num_dry = len(rain_df.index)-num_wet
    else:
        num_wet = ind.sum(axis=0).values[0]
        num_dry = len(rain_df.index)-num_wet-1


    #takes the difference between day and previous day indicator then set T or F if equal to 1
    wd = ind.diff(axis=0).eq(1)
    dw = ind.diff(axis=0).eq(-1)

    #if diff equals to 1, then that means that day is a wet day w/ previous day being a dry day
    num_wd = wd.apply(lambda x: [1 if y==True else 0 for y in x])
    num_dw = dw.apply(lambda x: [1 if y==True else 0 for y in x])
    #sum the number of wet dry day sequences and dry wet day sequences
    num_wd = num_wd.sum(axis=0).values[0]
    num_dw = num_dw.sum(axis=0).values[0]
    #if there are no dry days or wet days
    if num_dry==0:
        pwd = None
    else:
        pwd = num_wd/num_dry

    if num_wet == 0:
        pdw = None
    else:
        pdw = num_dw/num_wet
    return pwd,pdw

#find average pwd and pdw of january
pdw_av=0
pwd_av=0

wet_yr = 0
dry_yr = 0
for yr in range(200):
    january_all = og_data[og_data['month']==1]
    january_one = january_all[january_all['year']==2002+yr]
    rain = january_one[['pred']]

    pwd, pdw = gen_monthly_probs(rain)

    if pwd is not None:
        dry_yr+=1
        pwd_av+=pwd

    if pdw is not None:
        wet_yr+=1
        pdw_av+=pdw

pdw_av=pdw_av/wet_yr
pwd_av=pwd_av/dry_yr

def circ(i,l,u):
    if i < l:
        i=i+(u-l+1)
    else:
        i=i-(u-l+1)
    return i

def linint(month,y,xp):
    yp = []
    cc=month-1
    prv=circ(cc,1,12)
    x = month-1.5
    difx = xp-x
    if difx > 12:
        difx -= 12
    for j in range(14):
        val = y[prv][j]+(y[month][j]-y[prv][j])*difx
        yp.append(val)
    return yp


def cfact(y):
    values = {'rain':0,'tnw':0,'tnd':0,'txw':0,'txd':0,'srad':0}
    cf = {1:values,2:values,3:values,4:values,5:values,6:values,7:values,8:values,9:values,10:values,11:values,12:values}
    for i in cf:
        for j in range(10):
            cc = j - 1
            prv = circ(cc,1,10)
            cc=i+1
            nxt = circ(cc,1,10)




#rnum of only january
rnum = input.loc[month-1,'rnum']
#new alpha and beta
#alpha = min(max(a,0.01),0.998)
alpha = a
beta = b
#beta = max(rnum/alpha,0)
print(np.random.gamma(alpha,beta))
#print(np.random.gamma(a,b))

print(alpha,beta)

