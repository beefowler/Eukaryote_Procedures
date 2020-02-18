%This is a script to test whether our volume bins are well resolved
%Prompted by the comments of Reviewer 2, we will test to see that the number of volume bins we selected correctly resolves the dynamics. 
    %We will do this by 
        %Choosing a random subset of days from different times of year
        %Doubling the Volbins vector so that we have twice the resolution 
        %Recreating the model inputs according to that Volbins
        %Fitting (a modified, half-timestep) model to the new input. 
        %Comparing the simulation dynamics and division rate estimates across the two models. 

        
%1 ) We will test on 2017, since that is our most recent complete year of data

%a = 736785 %first day of 2017 data
%b = 737059 %last day of 2017 data
%days = floor((b-a).*rand(10,1) + a)



% 2) Becuase the eukaryote volbins we used was 2.^[-5:1/5:8]. We will just
% switch 1/5 to 1/10. 

%double_volbins = 2.^[-5:1/10:8] 

%the new vector is 1x131 and spans the same range of sizes
%save('Double_volbins.mat', 'double_volbins', 'days')


% 3) Alright now time to use setup_days_picoeuks.m to generate the new
% higher-resolution inputs

%here are the inputs 
year2do = 2017;
beadpath = '///MVCO/FCB/MVCO_Jan2017/data/processed/beads/'; 
modelpath =   '///MVCO/FCB/MVCO_Jan2017/euk_model/doublebins_inputs/'; 
datapath = '//MVCO/FCB/MVCO_Jan2017/';
mergedpath0 = '//MVCO/FCB/MVCO_Jan2017/data/processed/grouped/merged/'; 
groupedpath =   '//MVCO/FCB/MVCO_Jan2017/data/processed/grouped/'; 
plotflag = 1;

%volbins = double_volbins; 

%now set up days 
%setup_days_picoeuks

%actually, this generates inputs for all the days of 2017, definitely
%overkill, but it was the simplest way to do this. 
%then went back in and removed all the unwanted days (not on days list) from the directory. 


% 4 ) Now we want to run ModelMVCO using our new inputs 

%filelist = dir('\\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\*data.mat'); 
%keyboard % double check that the above is the list of days you want to do! Not all the days. 
%filepath = '\\MVCO\FCB\MVCO_Jan2017\euk_model\doublebins_inputs\' ; 
%savepath = '\\MVCO\FCB\pico_euk_model\doublebin_outputs\';   

%ModelMVCO

% 5 ) now we see what happened. 
% 
outpath = '\\Overflow_Outputs_BLF\MVCO_Jan2017\'; 
new_outpath = '\\MVCO\FCB\pico_euk_model\doublebin_halftime_outputs3\'; 

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


%% Additional Analysis of Results  
% We'd like to generate figures that look at the two outputs 
compare_params = 1; %set to do or not

if compare_params == 1
for i = 1:length(days)
    
    cd 'C:\Users\blfow\Desktop\Eukaryote_Procedures\Before_Modeling'
    
    n = days(i)
    figure, clf 
    eval(['load ' modelpath 'day' num2str(n) 'data.mat'])
    eval(['load ' outpath 'day' num2str(n) 'output.mat'])
    
    load('Volbins.mat')
    load('Double_volbins.mat')

    theta = modelresults(2:15); 
    hr1 = 1; hr2 = 25; 
    N_dist = CONC; 

    [dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample(Einterp,N_dist,theta,volbins,hr1,hr2);
    dirsampledist = dirsample ./ sum(dirsample); 

   
    %make division and growth functions for subpopulation 1 
    helpful = volbins - volbins(5); %saves some typing
    del1=(theta(4).*helpful.^theta(2))./(1+(helpful.^theta(2))); 
    del1(1:5) = [0 0 0 0 0];
    y1=theta(1)*ones(size(Einterp)); 
    ind=find(Einterp < theta(3));
    y1(ind)=(theta(1)/theta(3)) * Einterp(ind); 

    %subpopulation 2 
    del2=(theta(8).*helpful.^theta(6))./(1+(helpful.^theta(6)));
    del2(1:5) = [0 0 0 0 0];
    y2=theta(5)*ones(size(Einterp));
    ind=find(Einterp < theta(7));
    y2(ind)=(theta(5)/theta(7)) * Einterp(ind);

    subplot(3,3,1)
        h = imagesc(1:25, 1:length(double_volbins), Vhists); 
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(Vhists)); 
        xlabel('Hour after dawn')
        ylabel('Cell size class')
        title(['Observed data'])

    subplot(3,3,2)
        h = imagesc(1:25, 1:length(volbins), simdist); 
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(simdist)); 
        xlabel('Hour after dawn')
        ylabel('Cell size class')
        title('Simulated distribution, original')


    subplot(3,3,4)
        plot(1:.5:length(volbins), 2*Vhists(:,hr1), 'color', [0.5 0.5 0.5])
        hold on 
        plot(1:length(volbins),Vt1(:,1),'color',[0 0.5 1], 'linewidth',1.5)
        plot(1:length(volbins),Vt2(:,1),'color',[0 0 0.8], 'linewidth', 1.5)
        title('Initial distributions')
        xlabel('Size class')
        ylabel('Proportion')  
        axis([0 66 0 .08])


    subplot(3,3,5)
        [~, iy] = sort(y1); 
        plot(Einterp(iy), y1(iy), '.-', 'color', [0 0.5 1], 'markersize', 8)
        [~,iy]=sort(y2);  set(gca,'box','on')
        hold on 
        plot(Einterp(iy),y2(iy),'.-','color',[0 0 0.8],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Fraction of cells growing')
        xlabel('Radiation (W m^{-2})')
        title('Growth functions')

    subplot(3,3,6)
        plot(1:length(volbins), del1, '.-', 'color', [0 0.5 1], 'markersize',8);
        hold on 
        plot(1:length(volbins),del2,'.-','color',[0 0 0.8],'markersize',8);

        if n == 736808
            axis([0 66, 0 0.001])
            text(10, 0.0006, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        elseif n < 736924
            axis([0 66, 0 0.02])
            text(7, 0.015, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        else 
            axis([0 66, 0 0.038])
            text(10, 0.028, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        end
        ylabel('Fraction of cells dividing')
        xlabel('Size class')
        title('Division functions')

        
    %now we'll do the outputs for the simulated data
        
    eval(['load ' new_outpath 'day' num2str(n) 'output.mat'])
       
    load('Volbins.mat')
    load('Double_volbins.mat')
    
    theta = modelresults(2:15); 
    hr1 = 1; hr2 = 25; 
    N_dist = CONC; 
    
    cd '\\\MVCO\FCB\pico_euk_model\EukWork\Half_timestep'
    [dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample(Einterp,N_dist,theta,double_volbins,hr1,hr2);
    dirsampledist = dirsample ./ sum(dirsample); 
    
    
    helpful = volbins - volbins(5);
    %make division and growth functions for subpopulation 1 
    del1=(theta(4).*helpful.^theta(2))./(1+(helpful.^theta(2))); 
    del1(1:5) = [0 0 0 0 0]; 
    y1=theta(1)*ones(size(Einterp)); 
    ind=find(Einterp < theta(3));
    y1(ind)=(theta(1)/theta(3)) * Einterp(ind); 

    %subpopulation 2 
    del2=(theta(8).*helpful.^theta(6))./(1+(helpful.^theta(6)));
    del2(1:5) = [0 0 0 0 0]; 
    y2=theta(5)*ones(size(Einterp));
    ind=find(Einterp < theta(7));
    y2(ind)=(theta(5)/theta(7)) * Einterp(ind);
    
            
    subplot(3,3,3)
        h = imagesc(1:25, 1:length(double_volbins), simdist); 
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(simdist)); 
        xlabel('Hour after dawn')
        ylabel('Cell size class')
        title('Simulated distribution, double bins') 
        
    subplot(3,3,7)
        plot(1:length(double_volbins), Vhists(:,hr1), 'color', [0 0.3 .6])
        hold on 
        plot(1:length(double_volbins),Vt1(:,1),'color',[1 0.5 0.5], 'linewidth',1.5)
        plot(1:length(double_volbins),Vt2(:,1),'color',[.6 .1 .2], 'linewidth', 1.5)
        title('Initial distributions')
        xlabel('Size class')
        ylabel('Proportion')    
        axis([0 131 0 .08])
        
    subplot(3,3,8)
        [~, iy] = sort(y1); 
        plot(Einterp(iy), y1(iy), '.-', 'color', [1 0.5 0.5], 'markersize', 8)
        [~,iy]=sort(y2);  set(gca,'box','on')
        hold on 
        plot(Einterp(iy),y2(iy),'.-','color',[.6 .1 .2],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Fraction of cells growing')
        xlabel('Radiation (W m^{-2})')
        title('Growth functions')
        if n == 736838
            axis([0 800 0 0.3])
        end

    subplot(3,3,9)
        plot(1:length(volbins), del1, '.-', 'color', [1 0.5 0.5], 'markersize',8);
        hold on 
        plot(1:length(volbins),del2,'.-','color',[.6 .1 .2],'markersize',8);

        if n == 736808
            axis([0 66, 0 0.001])
            text(10, 0.0006, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        elseif n < 736924
            axis([0 66, 0 0.02])
            text(7, 0.015, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        else 
            axis([0 66, 0 0.038])
            text(10, 0.028, {['\mu_1: ' num2str(round(modelresults(18)*100)/100)];...
                ['\mu_2: ' num2str(round(modelresults(19)*100)/100)];
                ['\mu : ' num2str(round(modelresults(17)*100)/100)]});
        end
        ylabel('Fraction of cells dividing')
        xlabel('Size class')
        title('Division functions')

        
    sgtitle(datestr(day))
    set(gcf, 'position', [8 2 936 800]) 
    
    savefig(['DoubleBinTest' datestr(day) '.fig'])
end
end
