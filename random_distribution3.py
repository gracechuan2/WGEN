import pandas as pd
import numpy as np
from scipy import stats
import math
import itertools
import time

#------------------------------------------------------------------------------------------------------------------------

#OFFICIAL NOTATION
#solar radiation = srad
#min temp = tmin
#max temp = tmax
#rainfall = rain
#min temp on dry days = tnd
#min temp on wet days = tnw
#max temp on dry days = txd
#max temp on wet days = txw

#NOTES
#producing 200 years of daily weather data takes 184.32431200000002 seconds
#coefficients are calculated on yearly basis as in they change per year as opposed to using the same coefficients across
# all years in the input data

#------------------------------------------------------------------------------------------------------------------------

#MANUAL INPUTS
n_yrs= 1
#start year doesn't super matter to be honest, just to keep notation consistent and easy to compare
strt_yr = 2002

#------------------------------------------------------------------------------------------------------------------------

#function to create lag 1 serial correlation matrix
def lag_corr(data,lag):
    #compute the diagonal
    diagn = []
    for j in ['tmnd','tmxd','solar']:
        row = []
        for k in ['tmnd','tmxd','solar']:
            values1 = data.loc[j,:]
            values2 = data.loc[k,:]
            n=len(values1)
            y1 = values1[lag:]
            y2 = values2[:n-lag]
            corr = np.corrcoef(y1, y2, ddof=0)[0, 1]
            row.append(corr)
        diagn.append(row)
    return np.array(diagn)

def transition_matrix(order,params):
    pwd = params['pwd']
    pdd = 1-pwd
    pdw = params['pdw']
    pww = 1-pdw
    #our two states, wet and dry
    states = ['w','d']
    queue = ['w','d']
    final_probs = {}
    #BFS to get all combos of states with length 3
    while queue!=[]:
        ele = queue.pop(0)
        if len(ele) ==order:
            queue.append(ele)
            break
        for s in states:
            new= ele+s
            queue.append(new)
    #create probability matrix
    for state in queue:
        p = 1
        p_list = []
        for i in range(len(state)-1):
            if state[i] == 'w' and state[i+1]=='w':
                p*=pww
            elif state[i]=='w' and state[i+1]=='d':
                p*=pdw
            elif state[i]=='d' and state[i+1]=='d':
                p*=pdd
            else:
                p*=pwd
        p1=p
        p2=p
        for i in states:
            if i=='w':
                if state[-1]=='w':
                    p1*=pww
                else:
                    p1*=pwd
            else:
                if state[-1]=='w':
                    p2*=pdw
                else:
                    p2*=pdd
        #reweighting probabilities
        total = p1+p2
        p1 = p1/total
        p2 = p2/total
        p_list.append(p1)
        p_list.append(p2)
        final_probs[state] = p_list

    return final_probs


class WGEN:
    #input_mns_sds and input_data are strings
    def __init__(self, input_mns_sds,input_data, strt_yr,n_yrs,xlat):
        self.n_yrs =n_yrs
        self.strt_yr = strt_yr
        self.means_sds = {}
        self.input = pd.read_excel(input_mns_sds)
        self.xlat = xlat
        #INPUT existing weather dataset
        self.og_data = pd.read_csv(input_data)
        #keep avg of generated monthly precipitation values and avg
        self.generated_avgs = {}
        #convert table of means and sds for each month into easily accessible dictionary:
        for month in range(12):
            month+=1
            #each month will have a weather dict with 8 keys for each weather variable (rain, min temp, max temp, srad)
            weather_dict = {'rain':[self.input.loc[month-1,'rain1']/self.input.loc[month-1,'rnum'],self.input.loc[month-1,'rain2']], 'tmin':[self.input.loc[month-1,'tmin'],self.input.loc[month-1,'tmin2']], 'tmax':[self.input.loc[month-1,'tmax'],self.input.loc[month-1,'tmax2']], 'srad':[self.input.loc[month-1,'srad'],self.input.loc[month-1,'srad2']], 'tnw':[self.input.loc[month-1,'tnmnw'],self.input.loc[month-1,'tnsdw']],'tnd':[self.input.loc[month-1,'tnmnd'],self.input.loc[month-1,'tnsdd']],'txw':[self.input.loc[month-1,'txmnw'],self.input.loc[month-1,'txsdw']],'txd':[self.input.loc[month-1,'txmnd'],self.input.loc[month-1,'txsdd']]}
            self.means_sds[month] = weather_dict
            self.generated_avgs[month] = {'rain':0,'tnw':0,'tnd':0,'txw':0,'txd':0,'srad':0}
        #M0 is the correlation matrix of variables on the same day, M1 is same as M0 but one variable is lagged
        self.weather = self.og_data[['tmnd','tmxd','solar']].T
        M0=np.corrcoef(self.weather)
        M1 = lag_corr(self.weather,1)
        #calculate A and B matrices
        self.A = M1@np.linalg.inv(M0)
        BB = M0-M1@np.linalg.inv(M0)@M1.T
        # SOLVE FOR B
        eigens = np.linalg.eigh(BB)
        P = eigens[1]
        D = np.diag(eigens[0])
        self.B = P@D**(1/2)

    def gen_monthly_probs(self,rain_df):

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

    #generating probabilities of wet day given dry and dry day given wet; averages the probabilities across years per month
    def gen_rain_probs(self,data):
        #generate probabilities of wet day given a wet day, wet day  given a  dry day, etc. for the year
        prob_wd = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0}
        prob_dw = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0}

        for month in range(12):
            month+=1
            m = data[data['month']==month]
            wet_yrs = 0
            dry_yrs = 0
            for yr in range(self.n_yrs):
                m_yr = m[m['year']==self.strt_yr+yr]
                #single month given a specific yr
                #isolate the rainfall column
                rain_df = m_yr[['pred']]
                pwd, pdw = self.gen_monthly_probs(rain_df)

                if pwd is not None:
                    dry_yrs+=1
                    prob_wd[month]+=pwd

                if pdw is not None:
                    wet_yrs+=1
                    prob_dw[month]+=pdw

            #divide sums
            if wet_yrs!=0:
                avg = prob_dw[month]/wet_yrs
                prob_dw[month] = avg
            if dry_yrs!=0:
                prob_wd[month] = prob_wd[month]/dry_yrs

        prob_wd_df = pd.DataFrame({'month':list(prob_wd.keys()),'p(w|d)':list(prob_wd.values())})
        prob_dw_df = pd.DataFrame({'month':list(prob_dw.keys()),'p(d|w)':list(prob_dw.values())})

        #merge probability dataframes to original dataframe
        probs = pd.merge(prob_wd_df,prob_dw_df, on='month')

        return probs

    def weather_generation(self):
        probs = self.gen_rain_probs(self.og_data)
        data = pd.merge(probs,self.og_data,on='month')
        final_output = pd.DataFrame()
        #keep track of the years where it did rain
        wet_yrs = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0}
        for yr in range(self.n_yrs):
            current_yr = self.strt_yr+yr
            #make dataset smaller so its easier to work with
            data2 = data[data['year']==current_yr]
            doy=0
            #weather generation
            for m in range(12):
                m+=1
                month= data2[data2['month']==m].reset_index().drop(columns='level_0')
                #calculate beta and alpha for rain generation
                b = self.means_sds[m]['rain'][1]**2/self.means_sds[m]['rain'][0]
                a = self.means_sds[m]['rain'][0]/b
                pwd = month.loc[0,'p(w|d)']
                pdw = month.loc[0,'p(d|w)']
                params = {'a':a,'b':b,'pwd':pwd,'pdw':pdw}
                gen_rainfall = []
                gen_tntxsrad=[]
                w_or_d = None
                for d in month.index:
                    doy+=1
                    if d ==0:
                        w_or_d = 0
                        #solve for the first X0
                        sds = np.array([self.means_sds[m]['tnd'][1],self.means_sds[m]['txd'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
                        means = np.array([self.means_sds[m]['tnd'][0],self.means_sds[m]['txd'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
                        #find t1 for each newly generated yr
                        t1 = month[['tmnd','tmxd','solar']]
                        t1 = np.array(t1.iloc[0,:]).reshape(3,1)
                        X_i = (t1-means)/sds
                    #generate all weather variables
                    value_tntxsrad = self.gen_srad_temps(w_or_d,m,X_i,doy)
                    gen_tntxsrad.append(value_tntxsrad)
                    val, w_or_d  = self.rain_generation(params,w_or_d)
                    gen_rainfall.append(val)

                #concatenate resulting generated data on monthly basis
                gen_data = pd.concat(gen_tntxsrad,axis=1)
                gen_data.columns = [x for x in range(len(month.index))]
                #store away monthly rainfall total
                m_sum = sum(gen_rainfall)
                self.generated_avgs[m]['rain']+=(1/self.n_yrs)*m_sum
                #turn gen rainfall list into pd dataframe
                gen_rainfall = pd.DataFrame(gen_rainfall).T
                gen_data = pd.concat([gen_data,gen_rainfall],axis=0).T
                gen_data.columns = ['tmin','tmax','srad','rain']
                gen_data['year'] = self.strt_yr+yr
                gen_data['month'] = m
                gen_data['day'] = gen_data.index+1
                #find average of min and max temp
                #separate into 2 dataframes dry and wet
                dry = gen_data[gen_data['rain']==0]
                wet = gen_data[gen_data['rain']>0]

                tmnd_avg = float(dry['tmin'].mean())
                #check if there are no wet days, then the average would be just 0
                if tmnd_avg==float(gen_data['tmin'].mean()):
                    tmnw_avg = 0
                else:
                    tmnw_avg = float(wet['tmin'].mean())
                    #only need to count wet yrs once out of min and max temp
                    wet_yrs[m]+=1

                tmxd_avg = float(dry['tmax'].mean())
                if tmxd_avg==float(gen_data['tmax'].mean()):
                    tmxw_avg = 0
                else:
                    tmxw_avg = float(wet['tmax'].mean())

                #add the averages as sums first
                self.generated_avgs[m]['tnd']+=tmnd_avg
                self.generated_avgs[m]['tnw']+=tmnw_avg
                self.generated_avgs[m]['txd']+=tmxd_avg
                self.generated_avgs[m]['txw']+=tmxw_avg
                self.generated_avgs[m]['srad']+=float(gen_data['srad'].mean())/self.n_yrs
                #append newly generated data to growing dataset
                final_output=pd.concat([final_output,gen_data],axis=0)

        #divide yrs and number of wet yrs!
        for m in range(12):
            m+=1
            self.generated_avgs[m]['tnd']=self.generated_avgs[m]['tnd']/self.n_yrs
            self.generated_avgs[m]['tnw']=self.generated_avgs[m]['tnw']/wet_yrs[m]
            self.generated_avgs[m]['txd']=self.generated_avgs[m]['txd']/self.n_yrs
            self.generated_avgs[m]['txw']=self.generated_avgs[m]['txw']/wet_yrs[m]

        print(self.generated_avgs)
        return final_output.reset_index().drop(columns='index')

    #daily output
    def rain_generation(self,params, w_or_d):
        #list out the 4 weather parameters
        a = params['a']
        b = params['b']
        pwd = params['pwd']
        pdw = params['pdw']
        pww = 1-pdw
        pdd = 1-pwd
        #MARKOV CHAIN HERE
        #right now we use a first order MC! 1 is wet day and 0 is dry day
        if w_or_d==0:
            new_state =np.random.choice([1,0],replace=True,p=[pwd,pdd])
        else:
            if pww==0 and pdw ==0:
                new_state=0
            else:
                new_state =np.random.choice([1,0],replace=True,p=[pww,pdw])
        if new_state==1:
            return np.random.gamma(a,b), new_state
        else:
            return 0,new_state


    def calculate_residuals(self,X_i):
        #first compute epsilon, random independent components
        mu=0
        sigma=1
        e= np.random.normal(mu, sigma,3).reshape(3,1)
        X_i = np.dot(self.A,X_i)+np.dot(self.B,e)
        return X_i

    def gen_srad_temps(self,w_or_d,m,X_i,doy):
        X_i=self.calculate_residuals(X_i)
        #discern sds and means (this is where you would adjust the moving means)
        if w_or_d==0:
            sds = np.array([self.means_sds[m]['tnd'][1],self.means_sds[m]['txd'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
            means = np.array([self.means_sds[m]['tnd'][0],self.means_sds[m]['txd'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
        else:
            sds = np.array([self.means_sds[m]['tnw'][1],self.means_sds[m]['txw'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
            means = np.array([self.means_sds[m]['tnw'][0],self.means_sds[m]['txw'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
        t_i = X_i*sds+means
        #check if tmin > tmax, if this happens, swap the values; line 836 of fortran code;  also check if they are at least
        tmin = t_i[0][0]
        tmax = t_i[1][0]
        if tmin > tmax:
            t_i[0] = tmax
            t_i[1] = tmin
        #to ensure min and max are at least 0.1 apart
        if round(abs(tmin-tmax),3) <=0.1:
            t_i[1]+=1

        #incorporate the extraterrestrial radiation component
        srad = t_i[-1][0]
        rc = self.rada(doy)*0.8
        srmin=0.2*rc
        srad=min(max(srad,srmin),rc)
        t_i[-1]=srad
        return pd.DataFrame(t_i)

    def corrections(self):
        #monthly corrections
        rain_correct_factors = {}
        tnd_correct_factors = {}
        tnw_correct_factors = {}
        txd_correct_factors = {}
        txw_correct_factors = {}
        srad_correct_factors = {}
        #generate the dataset
        gen_dataset = self.weather_generation().reset_index()
        #add a rainy indicator
        gen_dataset['rainy'] = np.where(gen_dataset['rain']>0.01,1,0)

        for month in range(12):
            month+=1
            ind = self.generated_avgs[month]['rain']
            if ind==0:
                rain_correct_factors[month] =[1]
            else:
                rain_correct_factors[month]=[self.input.loc[month-1,'rain1']/ind]
            #tmin
            diff_tnd = self.means_sds[month]['tnd'][0]-self.generated_avgs[month]['tnd']
            diff_tnw = self.means_sds[month]['tnw'][0]-self.generated_avgs[month]['tnw']
            #tmax
            diff_txd = self.means_sds[month]['txd'][0]-self.generated_avgs[month]['txd']
            diff_txw = self.means_sds[month]['txw'][0]-self.generated_avgs[month]['txw']
            #srad
            diff_srad = self.means_sds[month]['srad'][0]-self.generated_avgs[month]['srad']

            #store away into libraries
            tnd_correct_factors[month] = [diff_tnd]
            tnw_correct_factors[month] = [diff_tnw]
            txd_correct_factors[month] = [diff_txd]
            txw_correct_factors[month] = [diff_txw]
            srad_correct_factors[month] = [diff_srad]

        #turn dictionaries into dfs
        rain_correct_factors = pd.DataFrame(rain_correct_factors).T
        rain_correct_factors.columns = ['rain_correction']
        rain_correct_factors['month'] = rain_correct_factors.index
        tnd_correct_factors = pd.DataFrame(tnd_correct_factors).T
        tnd_correct_factors['rainy'] = 0
        tnw_correct_factors = pd.DataFrame(tnw_correct_factors).T
        tnw_correct_factors['rainy'] = 1
        txd_correct_factors = pd.DataFrame(txd_correct_factors).T
        txd_correct_factors['rainy'] = 0
        txw_correct_factors = pd.DataFrame(txw_correct_factors).T
        txw_correct_factors['rainy'] = 1
        srad_correct_factors = pd.DataFrame(srad_correct_factors).T
        srad_correct_factors.columns = ['srad_correction']
        srad_correct_factors['month'] = srad_correct_factors.index

        tn_correct_factors = pd.concat([tnd_correct_factors,tnw_correct_factors],axis=0)
        tx_correct_factors = pd.concat([txd_correct_factors,txw_correct_factors],axis=0)

        tn_correct_factors=tn_correct_factors.rename(columns={0:'tn_correction'})
        tx_correct_factors=tx_correct_factors.rename(columns={0:'tx_correction'})

        #create columns for month in tn and tx dfs
        tn_correct_factors['month']  = tn_correct_factors.index
        tx_correct_factors['month']  = tx_correct_factors.index

        #merge correction dfs into main one
        gen_dataset  = pd.merge(gen_dataset,rain_correct_factors,on='month')
        gen_dataset = pd.merge(gen_dataset,tn_correct_factors,on=['month','rainy'])
        gen_dataset = pd.merge(gen_dataset,tx_correct_factors,on=['month','rainy'])
        gen_dataset = pd.merge(gen_dataset,srad_correct_factors,on=['month'])

        #apply corrections
        gen_dataset['tmin'] = gen_dataset['tmin']+gen_dataset['tn_correction']
        gen_dataset['tmax'] = gen_dataset['tmax']+gen_dataset['tx_correction']
        gen_dataset['rain'] = gen_dataset['rain']*gen_dataset['rain_correction']
        gen_dataset['srad'] = gen_dataset['srad']+gen_dataset['srad_correction']
        #sort by  index
        gen_dataset=  gen_dataset.sort_values(by=['index'])
        #increase index by 1
        gen_dataset['index']+=1
        return gen_dataset

        #extraterrestrial solar radiation calculation
    def rada(self,doy):
        #initialization
        t=2*math.pi*(doy+10)/365
        c1=math.cos(t)
        rad=math.pi/180
        #assign variables
        dec=-23.45*c1
        ssin=math.sin(rad*dec)*math.sin(rad*self.xlat)
        ccos=math.cos(rad*dec)*math.cos(rad*self.xlat)
        soc=ssin/ccos
        soc=min(max(soc,-1),1)
        dayl=12+24*math.asin(soc)/math.pi
        dsinb=3600*(dayl*ssin+24/math.pi*ccos*(1-soc**2)**(1/2))
        sc=1368*(1+0.033*math.cos(2*math.pi*doy/365))
        sod=sc*dsinb
        output=sod/10**6
        return output


#TESTING-------------------------------------------------------------------------------------------------------------------
t0=time.clock()
weather_gen = WGEN('input.xls','neutral_200.csv', 2002,200,19)
output = weather_gen.corrections()
t1=time.clock()-t0

print(t1)
output.to_csv('test outputs/python_generated12.csv',index=False)

#check if any

#export = output.to_csv('test outputs/full7.csv',index=False)


#print(t1)

#print(np.random.gamma(3,20))
#output.to_csv('test outputs/full3.csv',index=False)

#output=weather_gen.corrections()
#

#output.to_csv('output1.csv',index=False)

# december = weather_gen.og_data[weather_gen.og_data['month']==12]
# december1=december[december['year']==2002]['pred']
# december2=december[december['year']==2003]['pred']
#
# output = weather_gen.gen_monthly_probs(pd.DataFrame(december1))

#print(weather_gen.means_sds)
#print(final_output)

#------------------------------------------------------------------------------------------------------------------------




#export final datasets
# final_output['index']=final_output.reset_index().index+1
# final_output.to_csv('output1.csv', index=False)

#test means and sds

#print(gen_data)
#print(month[['tmnd','tmxd','solar','pred']])

#print(stats.mannwhitneyu(gen_data['rain'],month['pred']))

#DEBUGGING--------------
#
# #M0 is the correlation matrix of variables on the same day, M1 is same as M0 but one variable is lagged
# month2 = og_data[og_data['month']==10]
# #month2 = month2[month2['year']==2008].reset_index()
#
#
# weather2 = month2[['tmnd','tmxd','solar']].T
#
#
#
# M02=np.corrcoef(weather2)
# M12 = lag_corr(weather2,1)
#
# #calculate A and B matrices
# A2 = M12@np.linalg.inv(M02)
# BB2 = M02-M12@np.linalg.inv(M02)@M12.T
#
# # SOLVE FOR B, COME BACK
# eigens2 = np.linalg.eigh(BB2)
# P2 = eigens2[1]
# D2 = np.diag(eigens2[0])
#
# B2 = P2@D2**(1/2)
#
#Mann-Whitney U Tests
#import outputted dataset

def mann_whitneyu(x,y):
    p05 = 0
    p01 = 0
    p10 = 0
    p20 = 0
    for i in range(100):
        p = stats.mannwhitneyu(x,y)[1]
        if p > 0.05:
            p05+=1

        if p > 0.01:
            p01+=1

        if p > 0.10:
            p10+=1

        if p > 0.20:
            p20+=1

    p05 = p05/100*100
    p01 = p01/100*100
    p10 = p10/100*100
    p20 = p20/100*100

    return p05, p01, p10, p20



