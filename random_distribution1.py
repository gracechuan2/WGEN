import pandas as pd
import numpy as np
from scipy import stats
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
    def __init__(self, input_mns_sds,input_data, strt_yr,n_yrs):
        self.n_yrs =n_yrs
        self.strt_yr = strt_yr
        self.means_sds = {}
        self.input = pd.read_excel(input_mns_sds)
        #INPUT existing weather dataset
        self.og_data = pd.read_csv(input_data)
        #convert table of means and sds for each month into easily accessible dictionary:
        for month in range(12):
            month+=1
            #each month will have a weather dict with 8 keys for each weather variable (rain, min temp, max temp, srad)
            weather_dict = {'rain':[self.input.loc[month-1,'rain1']/self.input.loc[month-1,'rnum'],self.input.loc[month-1,'rain2']], 'tmin':[self.input.loc[month-1,'tmin'],self.input.loc[month-1,'tmin2']], 'tmax':[self.input.loc[month-1,'tmax'],self.input.loc[month-1,'tmax2']], 'srad':[self.input.loc[month-1,'srad'],self.input.loc[month-1,'srad2']], 'tnw':[self.input.loc[month-1,'tnmnw'],self.input.loc[month-1,'tnsdw']],'tnd':[self.input.loc[month-1,'tnmnd'],self.input.loc[month-1,'tnsdd']],'txw':[self.input.loc[month-1,'txmnw'],self.input.loc[month-1,'txsdw']],'txd':[self.input.loc[month-1,'txmnd'],self.input.loc[month-1,'txsdd']]}
            self.means_sds[month] = weather_dict
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


    #generating probabilities of wet day given dry and dry day given wet; averages the probabilities across years per month
    def gen_rain_probs(self,data):
        #generate probabilities of wet day given a wet day, wet day  given a  dry day, etc. for the year
        prob_wd = {}
        prob_dw = {}

        for month in range(12):
            month+=1
            m = data[data['month']==month]
            #isolate the rainfall column
            rain_df = m[['pred']]
            #creates an indicator column for if it rained that day or not
            ind = rain_df.apply(lambda x: [1 if y >0.01 else 0 for y in x]).reset_index().drop(columns = ['index'])
            #find the number of wet and dry days
            last_day = ind['pred'].iloc[-1]
            if last_day==1:
                num_wet = ind.sum(axis=0).values[0]-1
                num_dry = len(m.index)-num_wet
            else:
                num_wet = ind.sum(axis=0).values[0]
                num_dry = len(m.index)-num_wet-1
            #takes the difference between day and previous day indicator then set T or F if equal to 1
            wd = ind.diff(axis=0).eq(1)
            dw = ind.diff(axis=0).eq(-1)

            #if diff equals to 1, then that means that day is a wet day w/ previous day being a dry day
            num_wd = wd.apply(lambda x: [1 if y==True else 0 for y in x])
            num_dw = dw.apply(lambda x: [1 if y==True else 0 for y in x])

            num_wd = num_wd.sum(axis=0).values[0]
            num_dw = num_dw.sum(axis=0).values[0]

            if num_wet == 0:
                prob_dw[month] = num_dw
            else:
                prob_dw[month] = num_dw/num_wet

            prob_wd[month] = num_wd/num_dry

        prob_wd_df = pd.DataFrame({'month':list(prob_wd.keys()),'p(w|d)':list(prob_wd.values())})
        prob_dw_df = pd.DataFrame({'month':list(prob_dw.keys()),'p(d|w)':list(prob_dw.values())})

        #merge probability dataframes to original dataframe
        probs = pd.merge(prob_wd_df,prob_dw_df, on='month')

        return probs

    def weather_generation(self):
        probs = self.gen_rain_probs(self.og_data)
        data = pd.merge(probs,self.og_data,on='month')

        final_output = pd.DataFrame()
        for yr in range(self.n_yrs):
            #make dataset smaller so its easier to work with
            data = data[data['year']==self.strt_yr]
             #weather generation
            for m in range(12):
                m+=1
                month= data[data['month']==m].reset_index()
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
                    if d ==0:
                        w_or_d = 0
                        #solve for the first X0
                        sds = np.array([self.means_sds[m]['tmin'][1],self.means_sds[m]['tmax'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
                        means = np.array([self.means_sds[m]['tmin'][0],self.means_sds[m]['tmax'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
                        t1 = np.array(self.weather.iloc[:,0]).reshape(3,1)
                        X_i = (t1-means)/sds
                    gen_tntxsrad.append(self.gen_srad_temps(w_or_d,m,X_i))
                    val, w_or_d  = self.rain_generation(params,w_or_d)
                    gen_rainfall.append(val)

                #concatenate resulting generated data on monthly basis
                gen_data = pd.concat(gen_tntxsrad,axis=1)
                gen_data.columns = [x for x in range(len(month.index))]
                gen_rainfall = pd.DataFrame(gen_rainfall).T
                gen_data = pd.concat([gen_data,gen_rainfall],axis=0).T

                gen_data.columns = ['tmin','tmax','srad','rain']
                gen_data['year'] = self.strt_yr
                gen_data['month'] = m
                gen_data['day'] = gen_data.index+1
                #append newly generated data to growing dataset
                final_output=pd.concat([final_output,gen_data],axis=0)

            self.strt_yr+=1

        return final_output

    #daily output
    def rain_generation(self,params, state):
        #list out the 4 weather parameters
        a = params['a']
        b = params['b']
        pwd = params['pwd']
        pdw = params['pdw']
        #construct the transition probability matrix
        matrix = transition_matrix(3,params)
        #MARKOV CHAIN HERE
        #COMEBACK -> what to do when its the first 3 days?
        if len(state)!=3:
            current_state = state[-1]
            if current_state=='d':
                p = pwd
            else:
                p = 1-pdw
            new_state =np.random.choice(['w','d'],replace=True,p=[p,1-p])
        else:
            new_state=np.random.choice(['w','d'],replace=True,p=matrix[state])
        #return new rain value from the gamma distribution and previous three states
        if new_state=='w':
            return np.random.gamma(a,b), state[1:]+'w'
        else:
            return 0,state[1:]+'d'

    def calculate_residuals(self,X_i):
        #first compute epsilon, random independent components
        mu=0
        sigma=1
        e= np.random.normal(mu, sigma,3).reshape(3,1)
        X_i = np.dot(self.A,X_i)+np.dot(self.B,e)
        return X_i

    def gen_srad_temps(self,w_or_d,m,X_i):
        X_i=self.calculate_residuals(X_i)
        #discern sds and means (this is where you would adjust the moving means)
        if w_or_d==0:
            sds = np.array([self.means_sds[m]['tnd'][1],self.means_sds[m]['txd'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
            means = np.array([self.means_sds[m]['tnd'][0],self.means_sds[m]['txd'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
        else:
            sds = np.array([self.means_sds[m]['tnw'][1],self.means_sds[m]['txw'][1],self.means_sds[m]['srad'][1]]).reshape(3,1)
            means = np.array([self.means_sds[m]['tnw'][0],self.means_sds[m]['txw'][0],self.means_sds[m]['srad'][0]]).reshape(3,1)
        t_i = X_i*sds+means

        return pd.DataFrame(t_i)

#weather_gen = WGEN('input.xls','neutral_200.csv', 2002,1)

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



