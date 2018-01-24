% Reads in statistics outputted by diablo 
turbo=dlmread('mean_res.txt');

% Set the domain size
NY=256;
LY=2*pi;
NYM=NY-1;

% Enter the viscosity
NU=0.01;
%PR=300;

% Enter the number of scalars
N_TH=2;

% Enter the richardson number for each scalar
RI(1)=64;

% Enter the nondimensional uptake rate
C0=0.0081;
% Enter the dimensional swimming speed
V_s=100e-6;
% Enter how often statistics were saved
save_stats=10;
C0=C0*mult^2;

% Nondimensional lengthscale
l_star=0.057/2/pi;

% Nondimensional Kolmogorov scale
eta=0.111;
% Dimensional Batchelor scale
l_b=eta*l_star/sqrt(PR);

% Set the starting time in code units for start of averaging
tstart=0;
tstart=0.3452;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1)))

row=1;
for k=1:nk
  tii(k)=turbo(row,1);
  thvar(k)=turbo(row,2);
  thme(k)=turbo(row,3);
  thme_0(k)=turbo(row,4);
  th_min(k)=turbo(row,5);
  th_max(k)=turbo(row,6);
  chi(k)=turbo(row,7);
  l_th(k)=turbo(row,8);
  flux(k)=turbo(row,9);
  row=row+1;
end

for k=1:nk
%l_test(k)=l_star*sqrt(thme(k)^2+thvar(k))/sqrt((chi(k)/2/NU)*PR);
l_test(k)=l_star*sqrt(thvar(k))/sqrt((chi(k)/2/NU)*PR);
l_c(k)=l_test(k);
end

for k=1:nk
%t_consumption(k)=thth(k*5,1)/(C0*(thth(k*5,2)+thth(k*5,1)))*mult^2;
t_mixing(k)=thvar(k)/chi(k)*mult^2;
t_chemotaxis(k)=l_c(k)/V_s;
t_uptake(k)=100;
end

l_c_mean=l_star*sqrt(mean(thvar+thme.^2)/mean(chi/2/NU*PR));
t_mixing_mean=mean(thvar)/mean(chi)*mult^2;
t_uptake_mean=100;
t_chemotaxis_mean=l_c_mean/V_s;
l_s=1e-3;
R1_mean=t_uptake_mean/t_chemotaxis_mean;
R2_mean=t_mixing_mean/t_chemotaxis_mean;
C_max=2.5e-3;
c_ratio=thvar(1)/(2*thme(1));

c_ratio*(l_s/l_c_mean)/(1/R1_mean+1/R2_mean)

