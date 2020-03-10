%In this script, we will analyze a subset of days with some typtical
%subpopulation behavior. 
%the goal will be to generate supplementary material that demonstrates the
%details of subpopulation simulations, parameters, and to
%use these days as examples that we refit the model to. 

outpath = '\MVCO_Jan2017\' ; %This is where we saved the original model outputs
inpath = '\MVCO_Jan2017\euk_model\dawnstart_inputs_2019\'; %this is where we saved the model inputs

plotting = 0; 
makesimulations = 0 ; 

%These are all from 2017. 
days =  [736808  % two distinct populations, don’t grow much
        736838  % two overlapping populations, one grows one doesn't really
        736873]; % one clear dist, one very wide, maybe to capture noise?

if plotting == 1
        
    %go through each day and plot
    for n = days'
    figure, clf 
    eval(['load ' inpath 'day' num2str(n) 'data.mat'])
    eval(['load ' outpath 'day' num2str(n) 'output.mat'])

    theta = modelresults(2:15); 
    hr1 = 1; hr2 = 25; 
    N_dist = CONC; 

    [~, simdist,Vt1,Vt2]=simdata_dirichlet_sample(Einterp,N_dist,theta,volbins,hr1,hr2);
    %get a sample from simdata as mnrnd rather than dirichlet. 
    % This elliminates a source of noise for our refitting process. 
    q = hr2-hr1; 
    dirsampledist=zeros(length(volbins),q);
    for i=1:q+1
        dirsampledist(:,i)=mnrnd(round(sum(N_dist(:,hr1-1+i))),simdist(:,i));
    end
    dirsampledist = dirsampledist ./ sum(dirsampledist);
        
    %make division and growth functions for subpopulation 1 
    helpful = volbins - volbins(5);
    del1=(theta(4).*helpful.^theta(2))./(1+(helpful.^theta(2))); 
    del1(1:5) = [0 0 0 0 0]; 
    y1=theta(1)*ones(size(Einterp)); 
    ind=find(Einterp < theta(3));
    y1(ind)=(theta(1)/theta(3)) * Einterp(ind); 

    %subpopulation 2 
    del2=(theta(8).*helpful.^theta(6))./(1+(helpful.^theta(6)));
    del2(1:5) = [0 0 0 0 0]; %our 2019 change which limits small cell division 
    y2=theta(5)*ones(size(Einterp));
    ind=find(Einterp < theta(7));
    y2(ind)=(theta(5)/theta(7)) * Einterp(ind);

    subplot(2,3,1)
        h = imagesc(1:25, 1:length(volbins), Vhists); 
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(Vhists)); 
        xlabel('Hour after dawn')
        ylabel('Cell volume (\mum^{3})')
        title(['Observed data'])

    subplot(2,3,2)
        h = imagesc(1:25, 1:length(volbins), simdist); 
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(simdist)); 
        xlabel('Hour after dawn')
        title('Simulated distribution')

    subplot(2,3,3)
        h = imagesc(1:25, 1:length(volbins), dirsampledist); 
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(dirsampledist)); 
        xlabel('Hour after dawn')
        title('Simulated data') 


    subplot(2,3,4)
        plot(1:length(volbins), Vhists(:,hr1), 'color', [0.5 0.5 0.5])
        hold on 
        plot(1:length(volbins),Vt1(:,1),'color',[0 0.5 1], 'linewidth',1.5)
        plot(1:length(volbins),Vt2(:,1),'color',[0 0 0.8], 'linewidth', 1.5)
        axis([0 66 0 .1])
        xticks([6 31 56])
        xticklabels([0.06 2 64])
        title('Initial distributions')
        xlabel('Cell volume (\mum^{3})')
        ylabel('Proportion')    

    subplot(2,3,5)
        [~, iy] = sort(y1); 
        plot(Einterp(iy), y1(iy), '.-', 'color', [0 0.5 1], 'markersize', 8)
        [~,iy]=sort(y2);  set(gca,'box','on')
        hold on 
        plot(Einterp(iy),y2(iy),'.-','color',[0 0 0.8],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Growth fraction')
        xlabel('Radiation (W m^{-2})')
        title('Growth functions')

    subplot(2,3,6)
        plot(1:length(volbins), del1, '.-', 'color', [0 0.5 1], 'markersize',8);
        hold on 
        plot(1:length(volbins),del2,'.-','color',[0 0 0.8],'markersize',8);
        xticks([6 31 56])
        xticklabels([0.06 2 64])
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

        ylabel('Division fraction')
        xlabel('Cell volume (\mum^{3})')
        title('Division functions')

    sgtitle(datestr(day))

    set(gcf, 'position', [93 89 952 529]) 
   savefig(['highres_' datestr(day) '.fig'])


    end %for n = days
end %if plotting == 1

%% Make Simulations
%Now we want to save the outputs of these simulations such that they can
% be used as inputs for model fitting
if makesimulations ==1
outpath = '\\MVCO_Jan2017\euk_model\outputs\' ; 
inpath = '\\MVCO_Jan2017\euk_model\dawnstart_inputs_2019\'; 
savepath = '\\pico_euk_model\Simulated_inputs_high_res\' ; %this is where our simulated inputs go 

for n = days'
%load inputs and relevant outputs
eval(['load ' inpath 'day' num2str(n) 'data.mat'])
load([outpath 'day' num2str(n) 'output.mat'], 'modelresults', 'CONC', 'Einterp')


%run simulation of best fit model 
theta = modelresults(2:15); 
hr1 = 1; hr2 = 25; 
N_dist = CONC; 

[~, simdist,Vt1,Vt2]=simdata_dirichlet_sample(Einterp,N_dist,theta,volbins,hr1,hr2);
    %let's try to just get a sample from simdata as mnrnd rather than
    %dirichlet 
    q = hr2-hr1; 
    dirsampledist=zeros(length(volbins),q);
    for i=1:q+1
        dirsampledist(:,i)=mnrnd(round(sum(N_dist(:,hr1-1+i))),simdist(:,i));
    end
    dirsampledist = dirsampledist ./ sum(dirsampledist);
    
%Saving simulations in data format 
Vhists = dirsampledist;    
N_dist = dirsampledist .* cellsperml; 

%save to Simulated_inputs/ directory 
save([savepath 'day' num2str(n) 'simd.mat'], 'day', 'volbins', 'Edata', 'Vhists', 'N_dist', 'cellsperml', 'dielstarthr') 

end 

end

%% Compare Results
%of data fit to observations, vs data fit to simulated observations based
%on best fit model

siminpath = '\\pico_euk_model\Simulated_inputs_high_res\'; 
simoutpath = '\\pico_euk_model\Simulated_outputs_high_res\'; 


%First we'll just compare the division rate estimates. 
dataframe = zeros(length(days), 2);

for i = 1:length(days)
    n = days(i); 
    eval(['load ' simoutpath 'day' num2str(n) 'output.mat'])
    
    subplot(length(days), 2, (i-1)*2+1)
    h = pcolor(CONC); 
    set(h, 'EdgeColor', 'none')
    title('Simulated Input') 
    
    subplot(length(days), 2,(i-1)*2+2) 
    h = pcolor(simPROPS); 
    set(h, 'EdgeColor', 'none');
    title('New Modelfit')
    
    dataframe(i, 2) = modelresults(17); 

    eval(['load ' outpath 'day' num2str(n) 'output.mat'])
    
    dataframe(i,1) = modelresults(17); 
    
end

figure
scatter(dataframe(:,1), dataframe(:,2))
hold on 
plot([0 2], [0 2])
ylabel('Estimated Division Rate (d^{-1})') 
xlabel('"True" Division Rate (d^{-1})')


%now let's look at parameter details for the new simulation 
%how about we just add a row to the previous figure. 

for i = 1:length(days)
    n = days(i); 
    figure, clf 
    eval(['load ' inpath 'day' num2str(n) 'data.mat'])
    eval(['load ' outpath 'day' num2str(n) 'output.mat'])

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
        h = imagesc(1:25, 1:length(volbins), Vhists); 
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        caxis([0 0.1])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(Vhists)); 
        xlabel('Hour after dawn')
        ylabel('Cell Volume (\mum^{3})')
        title(['Observed data'])

    subplot(3,3,2)
        h = imagesc(1:25, 1:length(volbins), simdist); 
        caxis([0 0.1])
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(simdist)); 
        xlabel('Hour after dawn')
        title('Simulated distribution')

    subplot(3,3,3)
        h = imagesc(1:25, 1:length(volbins), dirsampledist); 
        caxis([0 0.1])
        yticks([6 31 56])
        yticklabels([0.06 2 64])
        set(gca, 'Ydir', 'normal')
        set(h, 'AlphaData', ~isnan(dirsampledist)); 
        xlabel('Hour after dawn')
        title('Simulated data') 
        colorbar


    subplot(3,3,4)
        plot(1:length(volbins), Vhists(:,hr1), 'color', [0.5 0.5 0.5])
        hold on 
        plot(1:length(volbins),Vt1(:,1),'color',[0 0.5 1], 'linewidth',1.5)
        plot(1:length(volbins),Vt2(:,1),'color',[0 0 0.8], 'linewidth', 1.5)
        axis([0 66 0 .1])
        xticks([6 31 56])
        xticklabels([0.06 2 64])
        title('Initial distributions')
        xlabel('Cell volume (\mum^{3})')
        ylabel('Proportion')    

    subplot(3,3,5)
        [~, iy] = sort(y1); 
        plot(Einterp(iy), y1(iy), '.-', 'color', [0 0.5 1], 'markersize', 8)
        [~,iy]=sort(y2);  set(gca,'box','on')
        hold on 
        plot(Einterp(iy),y2(iy),'.-','color',[0 0 0.8],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Growth fraction')
        xlabel('Radiation (W m^{-2})')
        title('Growth functions')

    subplot(3,3,6)
        plot(1:length(volbins), del1, '.-', 'color', [0 0.5 1], 'markersize',8);
        hold on 
        plot(1:length(volbins),del2,'.-','color',[0 0 0.8],'markersize',8);
        xticks([6 31 56])
        xticklabels([0.06 2 64])
        
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
        ylabel('Division fraction')
        xlabel('Cell volume (\mum^{3})')
        title('Division functions')

        
    %now we'll do the outputs for the simulated data
        
    eval(['load ' siminpath 'day' num2str(n) 'simd.mat'])
    eval(['load ' simoutpath 'day' num2str(n) 'output.mat'])
       
    theta = modelresults(2:15); 
    hr1 = 1; hr2 = 25; 
    N_dist = CONC; 
    
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
        
    subplot(3,3,7)
        plot(1:length(volbins), Vhists(:,hr1), 'color', [0 0.3 .6])
        hold on 
        plot(1:length(volbins),Vt1(:,1),'color',[1 0.5 0.5], 'linewidth',1.5)
        plot(1:length(volbins),Vt2(:,1),'color',[.6 .1 .2], 'linewidth', 1.5)
        axis([0 66 0 .1])
        xticks([6 31 56])
        xticklabels([0.06 2 64])
        xlabel('Cell volume (\mum^{3})')
        ylabel('Proportion')    

    subplot(3,3,8)
        [~, iy] = sort(y1); 
        plot(Einterp(iy), y1(iy), '.-', 'color', [1 0.5 0.5], 'markersize', 8)
        [~,iy]=sort(y2);  set(gca,'box','on')
        hold on 
        plot(Einterp(iy),y2(iy),'.-','color',[.6 .1 .2],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Growth fraction')
        xlabel('Radiation (W m^{-2})')
        if n == 736838
            axis([0 800 0 0.3])
        end
           
        %these are for moving to side once made
         h = title({'Model Fit to'; 'Observed Data'}); 
         set(h,'Rotation',90)

    subplot(3,3,9)
        plot(1:length(volbins), del1, '.-', 'color', [1 0.5 0.5], 'markersize',8);
        hold on 
        plot(1:length(volbins),del2,'.-','color',[.6 .1 .2],'markersize',8);
        xticks([6 31 56])
        xticklabels([0.06 2 64])
        
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
        ylabel('Division fraction')
        xlabel('Cell volume (\mum^{3})')
        
        %these are for moving to side once made
        k = title({'Model Fit to'; 'Simulated Data'}) ; 
        set(k,'Rotation',90)
        
    sgtitle(datestr(day))
    set(gcf, 'position', [8 2 936 800]) 

    savefig(['refit_high_res' datestr(day) '.fig'])
end

