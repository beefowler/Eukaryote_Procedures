%Code used to calculate climatologies
% by assigning daily values to their corresponding day of the year 

filelist = dir('YearSummaries_Culled/*.mat');
filepath = 'YearSummaries_Culled/'; 

allconc = zeros(366, 16)*NaN; 
alldivrates = zeros(366, 16)*NaN; 
alllossrates = zeros(366, 16)*NaN; 
allnetmu = zeros(366, 16)*NaN; 

for i = 1:length(filelist)
    eval(['load ' filepath filelist(i).name])
    yearnum = year(datetime(divrate(1,1), 'convertfrom', 'datenum')); 
    jannum = datenum(['01-Jan-' num2str(yearnum)]); 
    y = yearnum - 2002; %y is essentially 'year of our time series' 1 being 2003
    
    datevals = divrate(:,1) - jannum + 1; %convert to day of year
    
    alldivrates(datevals,y) = divrate(:,2); 
    alllossrates(datevals, y) = lossrate(:,2);  
    allnetmu(datevals, y) = netmu(:,2); 
    
    allmatdatevals = allmatdate - jannum + 1;
    allmatdatevals = floor(allmatdatevals); 
    
 
    for d = 1:366;
    allconc(d, y) = nanmean(alleukconc(find(allmatdatevals==d))); %average concentration across each day first
    end   
    
    
end

%now average for each day of year across all the years in the time series
Climat_Conc = nanmean(allconc'); 
Climat_Div = nanmean(alldivrates'); 
Climat_Loss = nanmean(alllossrates'); 
Climat_Net = nanmean(allnetmu'); 

%save('Climatology_2019.mat', 'Climat_Div', 'Climat_Loss', 'Climat_Conc', 'Climat_Net')

