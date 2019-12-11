%Creates the projection matrix for a particular hour, given light data, 
%size classes, and relevant parameters.

%hr = hour we are focusing on, relevant if we want our model to
%restrict division before a certain time of day 
        %BE SURE TO CHANGE 'ts' TO REFLECT RESTRICTIONS ON MORNING DIVISION
        
%b = shape parameter of division function 
%E_star = shape parameter of growth function (point where function switches from linear to constant)
%Ymax = max fraction of cells growing into next size class

%volbins = designation of size classes of cells in population
%  eg. volbins = 2.^[x: DeltaV : y] 
%  DeltaV must be such that 1/DeltaV is an integer!
%Einterp = radiation data win W/m^2 



function [B] =matrix_const_plateau(hr,Einterp,volbins,b,dmax,E_star,Ymax)

%Construct matrix A(t) for each time step within an hour based on delta and
%gammas at each 10 minute time intervals
%then construct B(t) which is A(t)s multiplied for the given hour
%multiply B(t)*w(t) to get the projection to next hour dist

%%%%%%--------Initial Parameters---------------
m=length(volbins);
j = find(2*volbins(1) == volbins);  %first index where division is allowed, cells twice as big as smallest bin. 
%j = find(2*volbins(1) <= volbins, 1); 
dt=(1/6); %corresponds to 10 mins (units are hours)
ts=0; %hours after dawn at which we will allow division

%-----------Delta function---------------------
%delta=frac of cells that divide
vstar = volbins(j-1); %last cell size that cannot divide. 
del=(dmax.*(volbins - vstar).^b)./(1+(volbins-vstar).^b); %subtract vstar from volbins to shift delta function over. 
del(1:(j-1)) = 0; %Want to limit division so smallest cells cant divide. 


%-----------Gamma function---------------------
%gamma (y) is frac of cells that grow into next class, depends on incident radiation (E), which is dependent on time
%Edata (first col time, 2nd col Irradiance):

y=Ymax*ones(size(Einterp));
ind=find(Einterp < E_star);
y(ind)=(Ymax/E_star) * Einterp(ind);

%-----------%constructing A---------------%
%%
A=spalloc(m,m,(m+2*(m-1))); %allocate sparse matrix for A, (size of matrix (m x m) and how many nonzeros will be in A
stasis_ind=1:(m+1):m^2;  %indexing for matrix -goes down a column and then right to next column
growth_ind=2:(m+1):m^2;  %corresponds to growth assignments
div_ind=((j-1)*m)+1:(m+1):m^2;  %division assignments 

%saveAs = {}; 

for t=1:(1/dt)      
     if hr < ts %before division is allowed
     delta=zeros(1,m);
     else
     delta=del; %allow division 
     end 
     A(stasis_ind)=(1-delta)*(1-y(t+6*hr)); % stasis, the 6*hr part in the indexing is because eahc hour is broken up inot 10min segments for the irradiance spline
     A(m,m)=(1-delta(m));
                 
     A(growth_ind)=y(t+6*hr)*(1-delta(1,1:m-1)); %growth on sub-diagonal
        
     A(1,1:j-1)=A(1,1:j-1)+2*delta(1:j-1); %division on superdiagonal and on first row partial
     A(div_ind)=2*delta(j:m);
     %saveAs = [saveAs {A}];

     %multiply A's for each ten minute interval to get B for the hour 
     if t ==1
         B=A;
     else
        B=A*B;
     end
     A = A.*0;
     if isnan(B),
         disp('B has a NaN!')
%                 keyboard, 
     end
     
end


    