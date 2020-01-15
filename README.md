# Eukaryote_Procedures
A clean version of the tools used for estimating division rates of picoeukaryotes, including descriptions and instructions. 

Here are the tools used to generate the results in Fowler et al. (in review). 

The general procedures are as follows, and files in this repository are grouped into folders accordingly. 
  1. Before_Modeling- prepare data inputs for model. 
  2. Model_Essential_Scripts- optimize model to observed data. 
  3. Post_Modeling- synthesize results and filter. 
  4. Figures_Data- data packaged to generate the figures in paper. 


# Environment 
    MATLAB 
    MATLAB Optimization Toolbox 
    MATLAB Parallel Computing Toolbox

Our model optimizations were all run on MATLAB R2018a, but R2019a also seems to work just fine. 


# Before Modeling
Here are the steps we took to generate the model inputs from the raw FCB data. Please contact us if you would like to work from the raw data. 

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

This generates a .mat file for each day of the year for which we have data. In it are the variables: \ 
day - the matlab date number \ 
volbins - just useful for backtracking \ 
dielstarthr - the (rounded) number of hours after midnight at which dawn occurs  \ 
cellsperml - 1x25 with concentrations of eukaryotes at every hour begining at dawn \ 
Edata - 2 x 73 matrix. First column is # of hours after dawn, second column is incidant radiation (W m^{-2}) \ 
N_dist - 66 x 25 matrix. Counts of cells in each bin for each hour of the day, beginning at dawn. \ 
Vhists - 66 x 25 matrix. Proportion of cells in each bin for each hour of the day (each hour sums to 1). \ 


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
Einterp - Interpolated light data for that day, values in (W m^{-2}) every 10 minutes begining at dawn. \ 
CONC - 66 x 25 matrix. Observed concentrations of cells in each bin for each hour of the day, beginning at dawn. \ 
simCONC - 66 x 25 matrix. Simulated concentrations of cells in each bin for each hour of the day according to best fit parameters. \ 
simPROPS - 66 x 25 matrix. Simulated proportion of cells in each bin for each hour of the day according to best fit parameters. Each hour sums to 1. \ 
modelfits - X by 17 matrix of all attempted model runs in the optimization process. Columns 1:14 are the parameters (details below). Column 15 has the negative log likelihood for each run. Column 16 is the division rate, and Column 17 is the ExitFlag from createOptimProblem for each run. \ 
allstarts - X by 14 matrix of the initial parameter values for each model run \ 
modelresults - 1 by 23 vector with results of optimization process. The first entry is the day in Matlab datenum form. Entries 2:15 are best-fit model parameters. modelresults(16) is the negative log likelihood. **modelresults(17) is the estimated division rate for the assemblage**, while modelresults(18:19) are the estimated division rates for each of the two subpopulations. modelresults(20:21) are the relative proportions of the two subpopulations at the end of the simulated day (as opposed to the starting proportions which is one of the parameters), and modelresults(23) is the ExitFlag from createOptimProblem for that best run. \ 


Througout our modeling code the vector of parameters, theta, is ordered as follows: 
    
    gmax1=theta(1); %max fraction of cells growing into next size class, subpopn 1 \
    b1=theta(2);  %shape parameter for division function, subpopn 1 \
    E_star1=theta(3); %shape parameter of growth function (point where function switches from linear to constant), subpopn 1 \ 
    dmax1=theta(4); %max fraction of cells able to divide in a given size class, subpopn 1 \ 
    gmax2=theta(5); %max fraction of cells growing into next size class, subpopn 2 \ 
    b2=theta(6); %shape parameter for division function, subpopn 2 \ 
    E_star2=theta(7); %shape parameter of growth function (point where function switches from linear to constant), subpopn 2 \ 
    dmax2=theta(8); %max fraction of cells able to divide in a given size class, subpopn 2 \ 
    f=theta(9); %proportion parameter, specifies starting fraction of subpopn 1 \ 
    m1=theta(10); %mean volume for starting cell size distribution, subpopn 1 \ 
    m2=theta(11); %mean volume for starting cell size distribution, subpopn 2 \ 
    sigma1=theta(12); %variance parameter for starting cell size distributions for popn 1 \ 
    sigma2=theta(13); %variance parameter for starting cell size distributions for popn 2 \ 
    s=theta(14); %overdispersion parameter for the Dirichlet-multinomial distribution \ 

