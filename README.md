# Eukaryote_Procedures
A clean version of the tools used for estimating division rates of Eukaryotes, including descriptions and instructions. 

Here are the tools used to generate the results in Fowler et al. (in review). 

The general procedures are as follows, and files in this repository are grouped into folders accordingly. 
  1. Before_Modeling- prepare data inputs for model. 
  2. Model_Essential_Scripts- optimize model to observed data. 
  3. Post_Modeling- synthesize results and filter. 
  4. Figures_Scripts- generate the figures in paper. 


# Prerequisites 
    MATLAB 
    MATLAB Optimization Toolbox 
    MATLAB Parallel Computing Toolbox

Our model optimizations were all run on MATLAB R2018a, but R2019a also seems to work just fine. 


# Before Modeling
Here are the steps we took to generate the model inputs from the FCB data, which are available [HERE].

Our value for the volbins vector, which is used to group cells into the size bins, can be found in Before_Modeling/Volbins.mat. 

For each year, we defined the paths to the appropriate directories of FCB data. 
    
    year2do = 2017;
    beadpath = '//MVCO_Jan2017/data/processed/beads/'
    modelpath =   '//MVCO_Jan2017/euk_model/dawnstart_inputs_2019/'
    datapath = '//MVCO_Jan2017/'
    mergedpath0 = '//MVCO_Jan2017/data/processed/grouped/merged/'
    groupedpath =   '//MVCO_Jan2017/data/processed/grouped/'
    plotflag = 1;

and then ran `setup_days_picoeuks.m.`

This generates a .mat file for each day of the year for which we have data. In it are the variables: 
day - the matlab date number 
volbins - just useful for backtracking 
dielstarthr - the (rounded) number of hours after midnight at which dawn occurs  
cellsperml - 1x25 with concentrations of eukaryotes at every hour begining at dawn 
Edata - 2 x 73 matrix. First column is # of hours after dawn, second column is incidant radiation (W m^{-2}) 
N_dist - 66 x 25 matrix. Counts of cells in each bin for each hour of the day, beginning at dawn. 
Vhists - 66 x 25 matrix. Proportion of cells in each bin for each hour of the day (each hour sums to 1). 


# Optimizing Model 
Optimizations for MVCO data were run one year at a time using the wrapper script ModelMVCO.m which can be found in the Model_Essential_Scripts/ 

Again we manually defined the paths to the relevant inputs and outputs. 
Usually something along these lines: 
    
    filelist = dir('\\MVCO_Jan2017\euk_model\dawnstart_inputs_2019\*data.mat');
    filepath = '\\MVCO_Jan2017\euk_model\dawnstart_inputs_2019\';
    savepath = '\\Outputs\MVCO_Jan2017\'; 

