# WGEN
Weather generation scripts based off of the Richardson (1984) climate model WGEN. Generates daily weather data for rainfall, minimum temperature, maximum temperature, and solar radiation from input daily weather data.
 
## How to Use
Inputs: 
- Input weather dataset (csv): daily values - sample input data set provided in repository 
   - Necessary variables: 
     - **"month"** month (1-12) 
     - **"year"** year 
     - **"day"** day (1-28,30,31) 
     - **"pred"** rainfall values (mm)
     - **"tmnd"** minimum temperature (°C)
     - **"tmxd"** maximum temperature (°C)
     - **"solar"** solar radiation (MJ/m**2/d)
 - Monthly means and sds (xls): monthly averages of each weather variable that the generated dataset should target - sample input file provided in repository
 *Use separate script to create this file from the input data or input values manually depending on your purposes* 
   - **"month"** month (1-12)
   - **"rainmn"** mean daily rainfall 
   - **"rainsd"** standard deviation of daily rainfall 
   - **"tminmn"** mean minimum temperature 
   - **"tminsd"** standard deviation of minimum temperature 
   - **"tmaxmn"** mean minimum temperature 
   - **"tmaxsd"** standard deviation of minimum temperature
   - **"sradmn"** mean solar radiation 
   - **"sradsd"** standard deviation of solar radiation 
   - **"tnmnd"**	mean minimum temperature for only dry days
   - **"tnsdd"** standard deviation of minimum temperature for only dry days
   - **"tnmnw"**	mean minimum temperature for only wet days
   - **"tnsdw"** standard deviation of minimum temperature for only wet days
   - **"txmnd"** mean maximum temperature for only dry days 
   - **"txsdd"** standard deviation of maximum temperature for only dry days
   - **"txmnw"**	mean maximum temperature for only wet days 
   - **"txsdw"** standard deviation of maximum temperature for only wet days
   - **"xmnd"**	mean solar radiation for only dry days
   - **"xsdd"** standard deviation of solar radiation for only dry days
   - **"xmnw"** mean solar radiation for only wet days
   - **"xsdw"** standard deviation of solar radiation for only wet days
 
 ## Quick Start
 1. Set repository WGEN as directory 
 2. Open "wgen.py" 
 3. In section 2, fill in the variables listed 
   3. input_mns_sds = name of the excel file that holds the input values for means and sds (make sure to type in with quotation marks 'string')
   3. input_data = name of the csv file with input dataset (make sure to type in with quotation marks 'string')
   3. strt_yr = start year (type in as a number)
   3. n_yrs = number of years you want to generate 
 4. Run script
 5. Resulting output data will be located in the "outputs" folder of the repository 

  
   
