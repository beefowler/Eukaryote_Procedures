%% Summarize Results for a Year
% SEVEN places to change year name before you start!

%%make a savename and Variable for grouped results
savename = 'ModelOutputs_2006.mat'; 
varname = 'AllResults_2006'; 
AllResults_2006 = []; 

filelist = dir('\\sosiknas1\Backup\Overflow_Outputs_BLF\MVCO_May2006_2019\da*output.mat');
filepath = '\\sosiknas1\Backup\Overflow_Outputs_BLF\MVCO_May2006_2019\'; 

%set up for making video
Writerobj1 = VideoWriter('ModelOutputs_2006.avi');
open(Writerobj1); 

obs_conc = zeros(1,25); 
obs_hr_mu = zeros(1,24);
pred_hr_mu = obs_hr_mu; 
hr_loss = obs_hr_mu;
pred_conc = obs_conc; 

for i = 1:length(filelist);
    eval(['load ' filepath filelist(i).name])
    
    %for 2006 only
    allstarts = allstarts{:}' ; 
    modelfits = modelfits{:}'; 
    
    
    AllResults_2006 = [modelresults; AllResults_2006]; %Save model results to grouped results
    
    %load the input variables too, they are usefull 
    day = modelresults(1); 
    inputpath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_May2006\euk_model\dawnstart_inputs_2019\'; 
    eval(['load ' inputpath 'day' num2str(day) 'data.mat'])

    %Interpolate light data
        time=0:(1/6):25;
        nnind = find(~isnan(Edata(:,2)));
        Edata=Edata(nnind,:);
        [unqE, eind]=unique(Edata(:,1));
        Einterp = interp1(Edata(eind,1),Edata(eind,2),time);
        Einterp(isnan(Einterp)) = 0;
        
    %make the figure
    Modelfit_Frame
    F1=getframe(gcf);%and add this frame to video
    
    writeVideo(Writerobj1, F1);
    clf
    
    sumconc = sum(CONC); 
    sumsim = sum(simCONC); 
    obs_conc(i,:) = sumconc; 
    obs_hr_mu(i,:) = 24*log(sumconc(2:end)./sumconc(1:(end-1))); 
    pred_hr_mu(i,:) = 24*log(sumsim(2:end)./sumsim(1:(end-1))); 
    pred_conc(i,:) = sumsim; 
    hr_loss(i,:) = pred_hr_mu(i,:) - obs_hr_mu(i,:); 
        

    
end

close(Writerobj1) 

%Flip it around cuz its upside down 
AllResults_2006 = AllResults_2006(end:-1:1,:);

save(savename, varname, 'obs_conc', 'obs_hr_mu', 'pred_hr_mu', 'pred_conc', 'hr_loss'); 


