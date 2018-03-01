% Reads in statistics outputted by diablo 
turbo=dlmread([rundir 'mean.txt']);

% Set the domain size
NY=64;
LY=6.28;
NYM=NY-1;

% Enter the viscosity
NU=0.00133;

% Enter the number of scalars
N_TH=1;

% Enter the richardson number for each scalar
RI(1)=1;

% Enter the nondimensional uptake rate
C0=0.0081;
% Enter the dimensional swimming speed
V_s=100e-6;
% Enter how often statistics were saved
save_stats=5;
% Enter the number of patches per side
mult=1;
C0=C0*mult^2;
% For low uptake simulation:
%C0=C0/10;
% For high uptake simulation:
%C0=C0*10;

%Set the initial nutrient variance
thvar_start=0.11931;

% Enter dimensional parameters for uptake calculations
% These are all fixed
phi=5e-12*1^3; %old
%phi=9.8657e-6;  %new in mgC
%phi=8.2214e-10;  %new in molC
lp=1e-3;
sigma=lp/2; %old
%sigma=2.5e-3; %new Patch size in meters
Cd=phi/((2*pi*sigma^2)^(3/2)); % Dimensional maximum nutrient concentration
Bd=5e11;
expenditure=1.138888e-11*(V_s^2);
L=0.057; %Domain size in meters
% For mid uptake
uptake=0.01/Bd %old, Uptake/cell
%uptake=0.01*(L/2/pi)^3;  %Uptake in total carbon
% For low uptake
%uptake=0.001/Bd;
% For high uptake
%uptake=0.1/Bd;

% Set the starting time in code units for start of averaging
%tstart=2.7379;
tstart=0.3452;
%tstart=16.035;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/2)

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  urms(k)=sqrt(turbo(row,1));
  vrms(k)=sqrt(turbo(row,2));
  wrms(k)=sqrt(turbo(row,3));
  ume(k)=turbo(row,1);
  vme(k)=turbo(row,2);
  wme(k)=turbo(row,3);
  epsilon(k)=turbo(row,1);
  row=row+1;
end

for k=1:nk
  tke(k)=0.5*(urms(k)^2+vrms(k)^2+wrms(k)^2);
end

% Now read in the geophysical variables
if (N_TH>0)
  turbo_th=dlmread([rundir 'mean_th.txt']);

% Determine the number of records in the file based on its length
  nk_th=floor(length(turbo_th(:,1))/(N_TH+1))
else
  nk_th=0;
end

nk_th=min(nk,nk_th)

row=1;
for k=1:nk_th
  tii(k)=turbo_th(row,2);
  dt(k)=turbo_th(row,3);
  row=row+1;
  for n=1:N_TH
    thth(k,n)=turbo_th(row,1);
    thvar(k,n)=turbo_th(row,2);
    thrms(k,n)=sqrt(thvar(k,n));
    row=row+1;
  end
end


% Extend thth using an exponential fit
k=nk_th
% Get the constant for the fit
A1=log(thth(nk_th,1))+2*C0*tii(nk_th);
A2=log(thth(nk_th,2))+2*C0*tii(nk_th);
while ((thth(k,1)>0.001*thth(1,1))||(thth(k,2)>0.001*thth(1,2)))
  k=k+1;
  tii(k)=tii(k-1)+dt(nk_th)*save_stats;
  thth(k,1)=exp(-2*C0*tii(k)+A1);
  thth(k,2)=exp(-2*C0*tii(k)+A2);
end    

% Integrate the uptake rates
kint_max=1;
tmax=60;
for k=1:nk
  if  ((tii(k)-tii(1))*mult^2*0.8108108<tmax) 
    kint_max=k;
  end
end

kint_max=length(tii);
%kint_max=nk_th;

start_time=tstart;
%start_time=tii(1);
%start_time=29.8951;

disp(['Time int (<bc>_m): ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,thth(1:kint_max,2)))]);
disp(['Time int(<bc>_nm): ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,thth(1:kint_max,1)))]);
disp(['Time int(<bc>_m-<bc>_nm): ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,thth(1:kint_max,2))-trapz((tii(1:kint_max)-tii(1))*mult^2*0.8108108,thth(1:kint_max,1)))]);

for k=1:length(thth(:,1))
  gain(k)=(thth(k,2)-thth(k,1))*uptake*Cd;
end

disp(['Motile uptake: ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,thth(1:kint_max,2)*uptake*Cd))]);
disp(['Non-Motile uptake: ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,thth(1:kint_max,1)*uptake*Cd))]);
disp(['Integrated Gain: ' num2str(trapz((tii(1:kint_max)-start_time)*mult^2*0.8108108,(thth(1:kint_max,2)-thth(1:kint_max,1))*uptake*Cd))]);

%h=area((tii(:)-tii(1))*mult^2*0.8108108,gain(:)-expenditure);

%plot((tii(5:nk_th)-tstart)*mult^2*0.8108108,gain(5:nk_th),'LineWidth',2);
%hold on
%plot((tii(nk_th:end)-tii(1))*mult^2*0.8108108,gain(nk_th:end),'-','LineWidth',2);
%h=line([0 (tii(end)-tii(1))*mult^2*0.8108108],[expenditure expenditure]);
%set(h,'LineStyle',':','LineWidth',2);
%set(gca,'FontName','Times');
%set(gca,'FontSize',14);
%xlabel('Time (sec)');
%title('Expenditure (dotted), Gain rate (solid, dashed) [mol_C/s]');
%title('Chemotactic Uptake Advantage [mol_C/s]');

%figure
%plot((tii(1:nk_th)-tstart)*mult^2*0.8108108,thth(1:nk_th,1)*uptake*Cd,'b-','LineWidth',2);
%plot((tii(1:nk_th)-tstart)*mult^2*0.8108108,(thth(1:nk_th,2)-thth(1:nk_th,1))*uptake*Cd,'b-','LineWidth',2);
%h=area((tii(1:end)-tstart)*mult^2*0.8108108/60,[thth(1:end,1)*uptake*Cd (thth(1:end,2)-thth(1:end,1))*uptake*Cd]);
%plot((tii(1:nk)-tstart)*mult^2*0.8108108/60,thth(1:nk,1)*uptake*Cd,'b-','LineWidth',2);
%hold on
%plot((tii(1:nk)-tstart)*mult^2*0.8108108/60,thth(1:nk,2)*uptake*Cd,'r-','LineWidth',2);
%plot((tii(nk:kint_max)-tstart)*mult^2*0.8108108/60,thth(nk:kint_max,1)*uptake*Cd,'b--','LineWidth',2);
%plot((tii(nk:kint_max)-tstart)*mult^2*0.8108108/60,thth(nk:kint_max,2)*uptake*Cd,'r--','LineWidth',2);
%plot((tii(1:nk_th)-tstart)*mult^2*0.8108108,thth(1:nk_th,2)-thth(1:nk_th,1),'r-')

%plot((tii(nk_th:end)-tii(1))*mult^2*0.8108108,thth(nk_th:end,2)*uptake*Cd,'r--','LineWidth',2);
%h=line([0 (tii(end)-tii(1))*mult^2*0.8108108],[expenditure expenditure]);
%set(h,'LineStyle',':','LineWidth',2);
%set(gca,'FontName','Times');
%set(gca,'FontSize',14);
%xlabel('Time (sec)');
%title('Expenditure (dotted), Uptake (red=motile, blue=non-motile) [mol_C/s]');



