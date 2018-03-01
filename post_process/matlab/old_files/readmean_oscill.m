% This version of readmean.m is designed for oscillating flows
% Reads in statistics outputted by diablo 
turbo=dlmread('mean.txt');

% Set the domain size
NY=32;
LY=1;
NYM=NY-1;

% Enter the viscosity
NU=1/400;

% Enter the number of scalars
N_TH=0;

% Enter the richardson number for each scalar
RI(1)=100;

% Set the starting time in code units for start of averaging
tstart=8;

% The number of bins to user per cycle
ncycle=12;
% The frequency of oscillation
omega0=0.314;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/(NY+3))


row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  ubulk(k)=turbo(row,1);
  row=row+1;
  for j=1:NY
    gyf(j)=turbo(row,2);
    ume(j,k)=turbo(row,3);
    vme(j,k)=turbo(row,4);
    wme(j,k)=turbo(row,5);
    urms(j,k)=turbo(row,6);
    vrms(j,k)=turbo(row,7);
    wrms(j,k)=turbo(row,8);
    uv(j,k)=turbo(row,9);
    uw(j,k)=turbo(row,10);
    wv(j,k)=turbo(row,11);
    dudy(j,k)=turbo(row,12);
    dwdy(j,k)=turbo(row,13);
    cp(j,k)=turbo(row,14);
    shear(j,k)=turbo(row,15);
    row=row+1;
  end
end


% Now read in the geophysical variables
if (N_TH>0)
  turbo_th=dlmread('mean_th.txt');

% Determine the number of records in the file based on its length
  nk_th=floor(length(turbo_th(:,1))/(N_TH*(NY+3)))
else
  nk_th=0;
end
row=1;
for k=1:nk_th
  tii(k)=turbo_th(row,2);
  dt(k)=turbo_th(row,3);
  row=row+1;
  ubulk(k)=turbo_th(row,1);
  row=row+1;
  for n=1:N_TH
  for j=1:NY
    thme(j,k,n)=turbo_th(row,3);
    dthdy(j,k,n)=turbo_th(row,4); % Add one here if a background was subtracted
    thrms(j,k,n)=turbo_th(row,5);
    thv(j,k,n)=turbo_th(row,6); 
    if (RI(n) ~= 0) 
      pe_diss(j,k,n)=turbo_th(row,7)*NU/RI(n);
    else
      pe_diss(j,k,n)=0;
    end
    row=row+1;
  end
  end
end



% Compute secondary quantities
for k=1:nk
  for j=1:NY
    tke(j,k)=0.5*(urms(j,k)^2.+vrms(j,k)^2.+wrms(j,k)^2.);
    if (dudy(j,k)~=0)
      nu_t(j,k)=-uv(j,k)/dudy(j,k);
    else
      nu_t(j,k)=0;
    end
% Calculate the vertical taylor scale
    if (shear(j,k)~=0) 
      taylor(j,k)=sqrt((ume(j,k)^2.+wme(j,k)^2.+urms(j,k)^2.+wrms(j,k)^2.)/shear(j,k));
    else
      taylor(j,k)=0;
    end
    if (N_TH > 0)
      for n=1:N_TH
        brunt(j,k,n)=sqrt(RI(n)*(dthdy(j,k,n))); 
        if (shear(j,k)~=0) 
          grarich(j,k,n)=brunt(j,k,n)^2./shear(j,k); 
        else
          grarich(j,k,n)=0;
        end
        if (dthdy(j,k,n)~=0)
          kappa_t(j,k,n)=-thv(j,k,n)/dthdy(j,k,n);
        else
          kappa_t(j,k,n)=0;
        end 
      end

    end
  end
end

% Get the time index based on start time
kstart=0;
for k=1:nk
  if (tii(k) <= tstart)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
'Start of time average: ',tii(kstart)

% Get an index for the time in terms of # of periods passed
period(1:nk)=tii(1:nk)*omega0/(2*pi);
count(1:ncycle)=0;
ume_mean(1:NY,1:ncycle)=0;
vme_mean(1:NY,1:ncycle)=0;
wme_mean(1:NY,1:ncycle)=0;
urms_mean(1:NY,1:ncycle)=0;
vrms_mean(1:NY,1:ncycle)=0;
wrms_mean(1:NY,1:ncycle)=0;
dudy_mean(1:NY,1:ncycle)=0;
dwdy_mean(1:NY,1:ncycle)=0;
tke_mean(1:NY,1:ncycle)=0;
uv_mean(1:NY,1:ncycle)=0;
wv_mean(1:NY,1:ncycle)=0;
cp_mean(1:NY,1:ncycle)=0;

for k=kstart:nk
  kindex=floor((period(k)-floor(period(k)))*ncycle)+1;
  count(kindex)=count(kindex)+1;
for j=1:NY
  ume_mean(j,kindex)=ume_mean(j,kindex)+ume(j,k);
  vme_mean(j,kindex)=vme_mean(j,kindex)+vme(j,k);
  wme_mean(j,kindex)=wme_mean(j,kindex)+wme(j,k);
  urms_mean(j,kindex)=urms_mean(j,kindex)+urms(j,k);
  vrms_mean(j,kindex)=vrms_mean(j,kindex)+vrms(j,k);
  wrms_mean(j,kindex)=wrms_mean(j,kindex)+wrms(j,k);
  dudy_mean(j,kindex)=dudy_mean(j,kindex)+dudy(j,k);
  dwdy_mean(j,kindex)=dwdy_mean(j,kindex)+dwdy(j,k);
  tke_mean(j,kindex)=tke_mean(j,kindex)+tke(j,k);
  uv_mean(j,kindex)=uv_mean(j,kindex)+uv(j,k);
  wv_mean(j,kindex)=wv_mean(j,kindex)+wv(j,k);
  cp_mean(j,kindex)=cp_mean(j,kindex)+cp(j,k);
  for n=1:N_TH
    thv_mean(j,n)=mean(thv(j,kstart:nk_th,n));
    dthdy_mean(j,n)=mean(dthdy(j,kstart:nk_th,n));
    thrms_mean(j,n)=mean(thrms(j,kstart:nk_th,n));
    thbar(j,n)=mean(thme(j,kstart:nk_th,n));
    pe_diss_mean(j,n)=mean(pe_diss(j,kstart:nk_th,n));
    if (dthdy_mean(j,n)~=0) 
      kappa_t_mean(j,n)=-thv_mean(j,n)/dthdy_mean(j,n);
    else
      kappa_t_mean(j,n)=0;
    end 
  end
  shear_mean(j)=mean(shear(j,kstart:nk));
end
end
for i=1:ncycle
  ume_mean(:,i)=ume_mean(:,i)/count(i);
  vme_mean(:,i)=vme_mean(:,i)/count(i);
  wme_mean(:,i)=wme_mean(:,i)/count(i);
end

  
%for j=1:NYM
%for k=1:nk
%  eta(j,k)=epsilon(j,k)^(-1/4)*(1/sqrt(400));
%end
%end
   
for j=2:NY
  gy(j)=(gyf(j)+gyf(j-1))/2;    
end
for j=2:NYM
  dyf(j)=(gy(j+1)-gy(j));
end
for j=2:NY
  dy(j)=gyf(j)-gyf(j-1);
end

for j=2:NY-1
  ry(j)=(gyf(j+1)-gyf(j))/(gyf(j)-gyf(j-1));
end


