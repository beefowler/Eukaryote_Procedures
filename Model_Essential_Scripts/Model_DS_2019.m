%This is for fitting the model to the Dilution Series (DS) data 

filelist = dir('DS_Analysis_2019\Inputs\*Data.mat');
filepath = 'DS_Analysis_2019\Inputs\'; 
savepath = 'DS_Analysis_2019\Outputs\'; 

%start while loop 
i = 1; 
while i <=length(filelist)
     
    %load days data 
    eval(['load ' filepath filelist(i).name])
    savename = [filelist(i).name(1:9) 'output.mat']; 
    
    N_dist = DS_PROPS; 
    Vhists = DS_COUNTS; 
    CONC = DS_CONC; 
    
    hr2 = size(N_dist, 2); 
        
        %apply model
        [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayMultiStart_hr2(day, Einterp, volbins, CONC, hr2);
        
        
        %save results
        save([savepath savename], 'CONC', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp', 'hr2')

    
   
   clearvars('-except', 'i', 'filelist', 'filepath', 'savepath'); 
      
i = i + 1; 
end



     