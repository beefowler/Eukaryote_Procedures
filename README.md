# Eukaryote_Procedures
This repository contains the tools, including corresponding descriptions and instructions, used for estimating the in situ division rates of picoeukaryotes at the Martha's Vineyard Coastal Observatory. 

The general procedures are as follows, and files in this repository are grouped into folders accordingly. 
  1. Before_Modeling- prepare data inputs for model. 
  2. Model_Essential_Scripts- optimize model to observed data. 
  3. Post_Modeling- trim and synthesize results. 
  4. Data_Products- data packaged to be analyzed or visualized, as in our 4 figures.  


# Environment 
    MATLAB 
    MATLAB Optimization Toolbox 
    MATLAB Parallel Computing Toolbox

Our model optimizations were all run on MATLAB R2018a, but R2019a also seems to work smoothly. 


# Before Modeling
In this section we describe the steps we took to generate the model inputs from the FCB data. The scripts used are available in the Before_Modeling directory. Please contact us if you would like to work from the model inputs or other upstream data products. 

Our value for the volbins vector, which is used to group cells into the size bins, can be found in Before_Modeling/Volbins.mat. If volbins is changed significantly, other modifications to the model scripts need to be considered. See below section on Model Sensitivity to Bin Size for more details. 

For each year, we ran setup_days_picoeuks.m, first defining the paths to the appropriate directories of FCB data. For example: 
    
    year2do = 2017;
    beadpath = '//MVCO_Jan2017/data/processed/beads/'
    modelpath =   '//MVCO_Jan2017/euk_model/dawnstart_inputs_2019/'
    datapath = '//MVCO_Jan2017/'
    mergedpath0 = '//MVCO_Jan2017/data/processed/grouped/merged/'
    groupedpath =   '//MVCO_Jan2017/data/processed/grouped/'
    plotflag = 1;
    
    setup_days_picoeuks

This generates a .mat file for each day of the year for which we have data. In it are the variables:   
day - the matlab date number   
volbins - the edges of our volume bins, here it is 1x66 vector.   
dielstarthr - the (rounded) number of hours after midnight at which dawn occurs on a day   
cellsperml - 1x25 vector with concentrations of eukaryotes at every hour beginning at dawn (cells/mL)   
Edata - 2 x 73 matrix. First column is # of hours after dawn, second column is incidant radiation (W m^{-2})   
N_dist - 66 x 25 matrix. Counts of cells observed in each volume bin for each hour of the day, beginning at dawn.   
Vhists - 66 x 25 matrix. Proportion of cells in each volume bin for each hour of the day (each hour sums to 1).   


Model Sensitivity to Bin Size     
Also in Before_Modeling/ is a folder called Test_BinResolution. Here, you can find the code that we used to test the sensitivity of our model to a change in the volume bin sizes. For a random subset of ten days in 2017, we redid the model setup described above with twice the number of volume bins over the same range of cell size. In order to properly compare this to the original model results, we needed to run our model with half the timestep (since cells can move up at most one size class per timestep, they would otherwise be prevented from growing as quickly as before). The code for that slighly modified model is in the Half_Timestep_Model folder, the test procedure is in the test_bin_resolution.m file and the results of our test can be seen in the ten .fig files. 

# Optimizing Model 
Optimizations for MVCO data were run one year at a time using the wrapper script ModelMVCO.m which can be found in the Model_Essential_Scripts/ folder. 

Again we manually defined the paths to the relevant inputs and outputs. 
Usually something along these lines: 
    
    filelist = dir('\\MVCO_Jan2017\euk_model\inputs_2019\*data.mat');
    filepath = '\\MVCO_Jan2017\euk_model\inputs_2019\';
    savepath = '\\Outputs\MVCO_Jan2017\'; 
    
For each day in the filelist, the script optimizes our model and creates an output with variabes as follows:   

  Einterp - Interpolated light data for that day, values in (W m^{-2}) every 10 minutes beginning at dawn.   
  CONC - 66 x 25 matrix. Observed concentrations of cells in each bin for each hour of the day, beginning at dawn.   
  simCONC - 66 x 25 matrix. Simulated concentrations of cells in each bin for each hour of the day according to best fit parameters.    
  simPROPS - 66 x 25 matrix. Simulated proportion of cells in each bin for each hour of the day according to best fit parameters. Each hour sums to 1.   
  modelfits - X by 17 matrix of all attempted model runs in the optimization process. Columns 1:14 are the parameters (details below). Column 15 has the negative log likelihood for each run. Column 16 is the division rate, and Column 17 is the ExitFlag from  createOptimProblem for each run.   
  allstarts - X by 14 matrix of the initial parameter values for each model run   
  modelresults - 1 by 23 vector with results of optimization process. The first entry is the day in Matlab datenum form. Entries 2:15 are best-fit model parameters. modelresults(16) is the negative log likelihood. **modelresults(17) is the estimated division rate for the assemblage**, while modelresults(18:19) are the estimated division rates for each of the two subpopulations. modelresults(20:21) are  the relative proportions of the two subpopulations at the end of the simulated day (as opposed to the starting proportions which is one of the parameters), and modelresults(22) is the ExitFlag from createOptimProblem for that best run. Modelresults(23) is the length of modelfits, or the number of model runs the optimization process tried.    


Throughout our modeling code the vector of parameters, theta, is ordered as follows:  
    
    gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1 
    b1=theta(2);  %shape parameter for division function, subpopn 1 
    E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1 
    dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1 
    gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2 
    b2=theta(6); %shape parameter for division function, subpopn 2 
    E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2 
    dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2 
    f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1 
    m1=theta(10); %mean volume for starting cell size distribution, subpopn 1 
    m2=theta(11); %mean volume for starting cell size distribution, subpopn 2  
    sigma1=theta(12); %variance parameter for starting cell size distributions for popn 1  
    sigma2=theta(13); %variance parameter for starting cell size distributions for popn 2  
    s=theta(14); %overdispersion parameter for the Dirichlet-multinomial distribution 

# After Modeling & Data Products
After the optimization process returned the best fit parameters and division rate estimates for each day, we organized the results by day and by year. For each year, we used the script GroupAYear2.m (which is available in the Post_Modeling directory) to generate a movie which compares the best-fit simulations to the observed data. Days for which the model was unable to reproduce the observed dynamics in the size distribution (1.1%; often due to apparent discontinuities in the size distribution) were trimmed. We used the script Climatology.m to group the results by day of the year in order to generate climatological averages. Our results are available in the Data directory. Data products are grouped by their corresponding figures in the manuscript and are described below. 

*All Model Results*     
The complete outputs from the Model Optimization process are available in the ModelResults_All.mat file. In this file, there is a matrix which contains the modelresults vector for each day of the time series with columns organized as described above. 

*Figure 1*  
Data_Fig1.mat includes the following variables:   
    alleukmatdate \[all eukaryotes matlab dates] = This is a 1x91722 vector of all the times of FCB observations in our time series. Times are in Matlab Datenum form UTC. Consecutive observations are roughly hourly.   
    alleukrunavg \[all eukaryotes running average] = This is the 1x91722 vector of concentration in cells/mL of picoeukaryotes observed at each time point, as a 48 hour running average. The values were calculated using the script, mvco_running_average.m which is available in the Model_Essential_Scripts directory, and which adjusts for large gaps in time in the dataset.   
    daily_divrates = This 3200x2 matrix contains the division rates for each day according to the results of our model-fitting process. The first column contains the days which had sufficient data for the model to be applied, again in Matlab Datenum form and the second column contains the estimated division rate (/day) for the picoeukaryote assemblage on each of those days. Fig1B, for example, would be generated by the command scatter(daily_divrates(:,1), daily_divrates(:,2)   
    daily_lossrates = this matrix is in the same form as daily_divrates, but the second column includes the calculated loss rates (division rate - observed net growth rate) for each day, rather than the division rates. The first column is the same as that of daily_divrates.  
To calculate climatologies, see information on Figure 2.    

*Figure 2*    
Note: Figure 2 describes climatologies for *Synechococcus* and Picoeukaryotes in terms of concentration, division rate, and estimates of Primary Productivity. Data_Fig2.mat includes the values for each day of the time series that are necessary for generating these climatologies. If only interested in the Climatological values, just take the averages of these matrices, e.g. plot(nanmean(daily_euk_conc')).     
The following variables are all matrices of size 366x16. The rows correspond to day of year and the columns are each distinct years, beginning with 2003 and ending with 2018. NaNs indicate days without data.          
    daily_euk_conc = average concentration of eukaryotes (cells/mL) observed on each day.     
    daily_euk_divrate = division rate (/day) estimate for picoeukaryote assemblage from model for each day.     
    daily_euk_lossrate = loss rate (/day) estimate for picoeukaryote assemblage for each day.     
    euk_min_cell_vol = the minimum value attained on each day by the mode of the picoeukaryote cell size distribution. To calculate primary productivity estimates, this value is assumed to be the "neutral" size which cells roughly grow from and divide back to.    
    
For each of the above, Data_Fig2.mat also includes the corresponding variables for *Synechococcus*, where "euk" is replaced with "syn". For *Synechococcus*, the minimum cell volume is taken from the dawn hour only, rather than looking for the minimum over the course of the day, because the mode of the Syn size distribution tends not to decrease after dawn.    

Lastly, for your convencience, we've also included the products of the Primary Productivity climatology estimates for picoeukaryotes (daily_euk_PP) and *Synechococcus* (daily_syn_PP) as 1x366 vectors containing values for each day of the year. These values are calculated from the above variables as described in the paper. 
    
*Figure 3*  
The data used to generate Figure 3 is in Data_Fig3.mat. As for Figure 2, we have chosen to make available all values for all days of the time series, though our figure only includes the averages for each year.   
All variables contain hourly rates in the form of a 16x366x24 array with the first dimension indicating the year (between 2003 and 2018), the second dimension is the day of the year, and the third is the hour of the day (beginning at dawn for that day).    
    alleuk_obs_hr_mu \[all eukaryotes, observed hourly net growth rate (mu)] = observed net growth rate (/day) in picoeukaryote assemblage calculated by comparing from each hour to the following hour.  
    alleuk_pred_hr_mu \[all eukaryotes, predicted hourly net growth rate] = hourly estimates of division rate (/day) calculated from the hourly changes in simulated population concentrations.   
    alleuk_hr_loss \[all eukaryotes, hourly loss rate] = hourly estimates of cell loss rate (/day) calculated by subtracting net growth rate from division rate.  
The corresponding variables for *Synechococcus* are available for each of the above, where "euk" is replaced with "syn". 

*Figure 4*    
The data used to generate Figure 4 in the manuscript. Once again, daily values are available for each day of the time series, while the figure only includes the climatological values. For each of the following matrices, each row is a year starting with 2003, and each column is the day of the year. NaNs indicate days without data.    
Daily_sunlight - 366x16 matrix with the mean incident radiation (W m^(-2)) for each day as measured by the PAR sensor at the MVCO Meteorological mast.    
Daily_temp - 366x16 matrix with the daily mean seawater temperature (&deg;C). Each row is a year starting with 2003 and each column is the day of the year. NaNs indicate days without data.     
daily_euk_divrate - 366x16 matrix with daily estimates of division rate for the picoeukaryote assemblage. Same as in Figure 2 data.    
More details about the environmental measurements are available at the [MVCO site](https://www.whoi.edu/website/mvco/meteorological-data/).     

Parameter Estimates from Known Values     
Also in the Post_Modeling folder is a script called subpop_case_studies_higher_res.m, which was used to generate Figures S9-S11 in the SI Appendix. This script focuses on three days, which were chosen to be "case studies" and are representative of the most common ways the two subpopulations are fit to the observed picoeukaryote dynamics at MVCO. The script plots the details of these subpopulations' behavior and parameter values, saves the optimized model outputs in the form of our model inputs, and then reruns the optimization process on these newly generated simulated data. The details of the second optimization are also plotted in the supplementary figures, such that we can compare the model output to the known population behavior and parameter values of a simulated data set. 


# Contact Us
If you have any questions about this repository or our manuscript, please email Bethany Fowler at bfowler@whoi.edu. 

