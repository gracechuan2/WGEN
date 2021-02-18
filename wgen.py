#Section 0---------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import math
import time

#Section 1------------------------------------------------------------------------------------------------------------------------
#Manual Input

input_mns_sds = 'input.xls'
input_data = 'neutral_200.csv'
strt_yr = 2002
n_yrs = 200
xlat = 19
corrections = True
output_name = 'full_gen4.csv'

#Section 2------------------------------------------------------------------------------------------------------------------------

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

class WGEN:
    #input_mns_sds and input_data are strings
    def __init__(self, input_mns_sds,input_data, strt_yr,n_yrs,xlat):
        self.n_yrs =n_yrs
        self.strt_yr = strt_yr
        self.means_sds = {}
        self.input = pd.read_excel(input_mns_sds)
        self.xlat = xlat
        #input existing weather dataset
        self.og_data = pd.read_csv(input_data)
        #keep avg of generated monthly precipitation values and avg
        self.generated_avgs = {}
        #convert table of means and sds for each month into easily accessible dictionary:
        for month in range(12):
            month+=1
            #each month will have a weather dict with 8 keys for each weather variable (rain, min temp, max temp, srad)
            weather_dict = {'rain':[self.input.loc[month-1,'rainmn'],self.input.loc[month-1,'rainsd']], 'tmin':[self.input.loc[month-1,'tminmn'],self.input.loc[month-1,'tminsd']], 'tmax':[self.input.loc[month-1,'tmaxmn'],self.input.loc[month-1,'tmaxsd']], 'srad':[self.input.loc[month-1,'sradmn'],self.input.loc[month-1,'sradsd']], 'tnw':[self.input.loc[month-1,'tnmnw'],self.input.loc[month-1,'tnsdw']],'tnd':[self.input.loc[month-1,'tnmnd'],self.input.loc[month-1,'tnsdd']],'txw':[self.input.loc[month-1,'txmnw'],self.input.loc[month-1,'txsdw']],'txd':[self.input.loc[month-1,'txmnd'],self.input.loc[month-1,'txsdd']]}
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
        self.prob_wd_df = None
        self.prob_dw_df = None
        #initialize the probability parameters through gen_rain_probs
        self.gen_rain_probs(self.og_data)
        #miscellaneous information needed
        self.num_days_months = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}

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
        self.prob_wd_df = prob_wd
        self.prob_dw_df = prob_dw

    def weather_generation(self):
        final_output = pd.DataFrame()
        #keep track of the years where it did rain
        wet_yrs = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0}
        #set t1 for Jan 1st of the start yr and find the very first X1
        t1 = np.array(self.weather.iloc[:,0]).reshape(3,1)
        sds = np.array([self.means_sds[1]['tnd'][1],self.means_sds[1]['txd'][1],self.means_sds[1]['srad'][1]]).reshape(3,1)
        means = np.array([self.means_sds[1]['tnd'][0],self.means_sds[1]['txd'][0],self.means_sds[1]['srad'][0]]).reshape(3,1)
        X_i = (t1-means)/sds
        w_or_d = 0
        for yr in range(self.n_yrs):
            doy=0
            #weather generation
            for m in range(12):
                m+=1
                #calculate beta and alpha for rain generation
                b = self.means_sds[m]['rain'][1]**2/self.means_sds[m]['rain'][0]
                a = self.means_sds[m]['rain'][0]/b
                pwd = self.prob_wd_df[m]
                pdw = self.prob_dw_df[m]
                params = {'a':a,'b':b,'pwd':pwd,'pdw':pdw}
                gen_rainfall = []
                gen_tntxsrad=[]
                for d in range(self.num_days_months[m]):
                    doy+=1
                    #generate all weather variables
                    X_i, value_tntxsrad = self.gen_srad_temps(w_or_d,m,X_i,doy)
                    gen_tntxsrad.append(value_tntxsrad)
                    val, w_or_d  = self.rain_generation(params,w_or_d)
                    gen_rainfall.append(val)
                #concatenate resulting generated data on monthly basis
                gen_data = pd.concat(gen_tntxsrad,axis=1)
                gen_data.columns = [x for x in range(self.num_days_months[m])]
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
                dry = gen_data[gen_data['rain']<=0.01]
                wet = gen_data[gen_data['rain']>0.01]
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
        return X_i, pd.DataFrame(t_i)

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
                rain_correct_factors[month]=[self.input.loc[month-1,'rainmn']*self.input.loc[month-1,'rnum']/ind]
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



#Section 3-------------------------------------------------------------------------------------------------------------------
#initialize weather generation class
weather_gen = WGEN('Inputs/'+input_mns_sds,'Inputs/'+input_data,strt_yr,n_yrs,xlat)
#run weather generation
if corrections:
    output = weather_gen.corrections()
else:
    output = weather_gen.weather_generation()
#export final dataset
output.to_csv('Outputs/'+output_name, index=False)
