%This is a script to test whether our volume bins are well resolved
%Prompted by the comments of Reviewer 2, we will test to see that the number of volume bins we selected correctly resolves the dynamics. 
    %We will do this by 
        %Choosing a random subset of days from different times of year
        %Doubling the Volbins vector so that we have twice the resolution 
        %Recreating the model inputs according to that Volbins
        %Fitting the model to the new input. 
        %Comparing the simulation dynamics and division rate estimates across the two models. 

        
%1 ) We will test on 2017, since that is our most recent complete year of data

%a = 736785 %first day of 2017 data
%b = 737059 %last day of 2017 data
%days = floor((b-a).*rand(10,1) + a)


% 2) Becuase there were 57 volume bins for Syn. let's make fewer eukaryote
% volbins. 

fewer_volbins = 2.^[-5:1/4:8]; 

%the new vector is 1x131 and spans the same range of sizes
save('fewer_volbins.mat', 'fewer_volbins', 'days')


% 3) Alright now time to use setup_days_picoeuks.m to generate the new
% lower-resolution inputs

%here are the inputs 
year2do = 2017;
beadpath = '//Sosiknas1/Lab_data/MVCO/FCB/MVCO_Jan2017/data/processed/beads/'; 
modelpath =   '//Sosiknas1/Lab_data/MVCO/FCB/MVCO_Jan2017/euk_model/fewerbins_inputs/'; 
datapath = '//Sosiknas1/Lab_data/MVCO/FCB/MVCO_Jan2017/';
mergedpath0 = '//Sosiknas1/Lab_data/MVCO/FCB/MVCO_Jan2017/data/processed/grouped/merged/'; 
groupedpath =   '//Sosiknas1/Lab_data/MVCO/FCB/MVCO_Jan2017/data/processed/grouped/'; 
plotflag = 1;

volbins = fewer_volbins; 

%now set up days 
setup_days_picoeuks

%actually, this generates inputs for all the days of 2017, definitely
%overkill, but it was the simplest way to do this. 
%then went back in and removed all the unwanted days (not on days list) from the directory. 


% 4 ) Now we want to run ModelMVCO using our new inputs 

%filelist = dir('\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\*data.mat'); 
%keyboard % double check that the above is the list of days you want to do! Not all the days. 
%filepath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\' ; 
%savepath = '\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\doublebin_outputs\';   

%ModelMVCO

% 5 ) now we see what happened. 

outpath = '\\sosiknas1\Backup\Overflow_Outputs_BLF\MVCO_Jan2017\'; 
new_outpath = '\\sosiknas1\Lab_data\MVCO\FCB\pico_euk_model\fewerbins_outputs\'; 

dataframe = zeros(10, 2); 

for i = 1:10
    n = days(i); 
    eval(['load ' outpath 'day' num2str(n) 'output.mat'])
    dataframe(i, 1) = modelresults(17); 
    
    subplot(10,3,3*(i-1)+1)
    h = pcolor(CONC); 
    set(h, 'EdgeColor', 'none')
    
    subplot(10,3,3*(i-1)+2)
    h = pcolor(simPROPS); 
    set(h, 'EdgeColor', 'none');
    
    eval(['load ' new_outpath 'day' num2str(n) 'output.mat'])
    dataframe(i, 2) = modelresults(17); 
    
    subplot(10,3,3*(i-1)+3)
    h = pcolor(simPROPS); 
    set(h, 'EdgeColor', 'none');
    
end

figure
scatter(dataframe(:,1), dataframe(:,2))
ylabel('Div Rate from Double Bins')
xlabel('Div Rate from Original Model')

