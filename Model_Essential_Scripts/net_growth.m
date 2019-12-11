function [n_mu] = net_growth(CONC, previousconc)

if exist('previousconc')
    n_mu = 24*log(sum(CONC(:,1)) / previousconc); 
else n_mu = NaN; 
end

for hour = 1:(size(CONC,2)-1)  
    n_mu = [n_mu 24*log(sum(CONC(:,hour+1)) / sum(CONC(:,hour)))]; 
end

end 
