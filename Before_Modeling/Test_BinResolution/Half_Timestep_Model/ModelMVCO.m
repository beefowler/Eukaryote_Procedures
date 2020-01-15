%% NOTE: THIS IS IN HALF_TIMESTEP FOLDER! Only use for double_bins test

filelist = dir('\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\*data.mat');
filepath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\';
savepath = '\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\doublebin_halftime_outputs3\'; 

%start while loop -
i = 1; 
while i <= length(filelist)
     
    %load days data 
    eval(['load ' filepath filelist(i).name])
    savename = [filelist(i).name(1:9) 'output.mat']; 
    
    
    if size(N_dist, 2) == 25
        
        %Interpolate light data
        time=0:(1/12):25;  %Interpolate to the correct time units
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
        
        %Change counts to concentrations
        CONC = Vhists .* cellsperml;
        
        %apply model
        [modelresults, modelfits, allstarts, simPROPS, simCONC] = OneDayMultiStart(day, Einterp, volbins, CONC);
        
        
        %save results
        save([savepath savename], 'CONC', 'modelresults', 'modelfits', 'allstarts', 'simPROPS', 'simCONC', 'Einterp')
               
    else
        disp([filelist(i).name(1:9) 'not enough hours'])
    end %if we have 25 hours
    
   
   clearvars('-except', 'Writerobj1', 'i', 'filelist', 'filepath', 'savepath'); 
   
i = i + 1; 

end

