#Section 0------------------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import xlwt

#Section 1------------------------------------------------------------------------------------------------------------------------
#Manual Input

input_filename = 'neutral_200.csv'

#Section 2------------------------------------------------------------------------------------------------------------------------
#set up table template and variables
table = {'month':[x+1 for x in range(12)]}

weather_vars = ['rainmn', 'tmaxmn', 'tminmn', 'sradmn', 'txmnd', 'txmnw', 'tnmnd', 'tnmnw', 'xmnd',	'xmnw',	'txsdd', 'txsdw', 'tnsdd', 'tnsdw',	'xsdd',	'xsdw',	'rainsd', 'tmaxsd',	'tminsd', 'sradsd','rnum']
for var in weather_vars:
    #make a new key for each weather variable
    table[var] = []

#import input data
input_data = pd.read_csv('Inputs/'+input_filename)
years = input_data['year'].unique()
#add rainfall indicator column
input_data['rainy']= np.where(input_data['pred']>0.01,1,0)

#iterate through each month and weather variable and compile
for m in range(12):
    m+=1
    #isolate dataframe
    month = input_data[input_data['month']==m]
    wet = month[month['pred']>0.01]
    dry = month[month['pred']<=0.01]
    #find means and sds of the 4 variables straightup without dry or wet days involved
    for weather in ['pred','solar','tmxd','tmnd']:
        #calculate mean and sds
        avg = month[weather].mean()
        sd = month[weather].std()
        #calculate the dry and wet day statistics too
        avg_dry = dry[weather].mean()
        sd_dry = dry[weather].std()
        avg_wet = wet[weather].mean()
        sd_wet = wet[weather].std()
        #add into dictionary
        if weather == 'pred':
            table['rainmn'].append(avg_wet)
            table['rainsd'].append(sd_wet)
        elif weather=='solar':
            table['sradmn'].append(avg)
            table['sradsd'].append(sd)
            table['xmnd'].append(avg_dry)
            table['xsdd'].append(sd_dry)
            table['xmnw'].append(avg_wet)
            table['xsdw'].append(sd_wet)
        elif weather=='tmxd':
            table['tmaxmn'].append(avg)
            table['tmaxsd'].append(sd)
            table['txmnd'].append(avg_dry)
            table['txsdd'].append(sd_dry)
            table['txmnw'].append(avg_wet)
            table['txsdw'].append(sd_wet)
        else:
            table['tminmn'].append(avg)
            table['tminsd'].append(sd)
            table['tnmnd'].append(avg_dry)
            table['tnsdd'].append(sd_dry)
            table['tnmnw'].append(avg_wet)
            table['tnsdw'].append(sd_wet)
    #calculate rnum
    rnum = 0
    for yr in years:
        month_yr = month[month['year']==yr]
        sum = month_yr['rainy'].sum()
        rnum+=month_yr['rainy'].sum()
    rnum=rnum/len(years)
    table['rnum'].append(rnum)

final_input = pd.DataFrame(table)
final_input = final_input.set_index('month')

#optional corrections of the means
def corrections(table):
    table_copy = table
    sds_vars = ['txsdd', 'txsdw', 'tnsdd', 'tnsdw',	'xsdd',	'xsdw',	'rainsd', 'tmaxsd',	'tminsd', 'sradsd']
    for weather in table.columns:
        if weather in sds_vars:
            for m in range(12):
                m+=1
                if m==1:
                    prv = 12
                    nxt = 2
                elif m==12:
                    nxt = 1
                    prv = 11
                else:
                    prv = m-1
                    nxt = m+1

                add = table_copy.loc[m,weather] - (0.25*(table_copy.loc[prv,weather]+2*table_copy.loc[m,weather]+table_copy.loc[nxt,weather]))

                #add correction to the weather variable
                if table_copy.loc[m,weather] + add > 0:
                    table.loc[m,weather] = table_copy.loc[m,weather]+add

    return table_copy

final_input = corrections(final_input)
#export final input table
final_input.to_excel('Inputs/input3.xls')





