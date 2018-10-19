% Read vertical profiles from mean.h5

% rundir='/local/scratch/public/cjh225/mini_diablo/Reb20aim/';
rundir='/store/DAMTP/cjh225/HPC_forcing_tests/256_Reb1/';

NU=6e-4; Re=1/NU; Ri_t=1; Pr=1;

fname=[rundir 'mean.h5'];

nk=h5readatt(fname,'/U1me','Samples');
NY=256;

U1me=zeros(nk,NY); U3me=zeros(nk,NY); THme=zeros(nk,NY);

for i=1:nk
    % Create the right dataset name for each time step
    if i<=10
        num=['000' int2str(i-1)];
    elseif i<=100
        num=['00' int2str(i-1)];
    elseif i<=1000
        num=['0' int2str(i-1)];
    else
        num=int2str(i-1);
    end

    % Read vertical profiles
    dname=['/U1me/' num];
    U1me(i,:)=h5read(fname,dname);
    dname=['/U3me/' num];
    U3me(i,:)=h5read(fname,dname);
    dname=['/THme/' num];
    THme(i,:)=h5read(fname,dname);
    dname=['/epsilon/' num];
    epsilon(i,:)=h5read(fname,dname);
    dname=['/chi/' num];
    chi(i,:)=h5read(fname,dname);
    dname=['/U1U2/' num];
    U1U2(i,:)=h5read(fname,dname);
    dname=['/U3U2/' num];
    U3U2(i,:)=h5read(fname,dname);
    dname=['/THflux/' num];
    THflux(i,:)=h5read(fname,dname);
end

% Compute kinetic energy in shear modes & associated viscous dissipation
EK_bar=0.5*mean(U1me.^2+U3me.^2,2);
CIKY = 1i*[0:ceil((NY-1)/2) ceil(-(NY-1)/2):-1];
CS1=fft(U1me,[],2);     CS1=CIKY.*CS1;      DU1DY=ifft(CS1,[],2,'symmetric');
CS1=fft(U3me,[],2);     CS1=CIKY.*CS1;      DU3DY=ifft(CS1,[],2,'symmetric');
S2=DU1DY.^2+DU3DY.^2;
eps_bar=1/Re*mean(DU1DY.^2+DU3DY.^2,2);

% Compute mean N2 & Richardson number
CS1=fft(THme,[],2);    CS1=CIKY.*CS1;      DTHDY=ifft(CS1,[],2,'symmetric');
N2=Ri_t*(1+DTHDY);
Ri_g=N2./S2;

% Compute buoyancy Reynolds number and shear production
Re_b=Re*epsilon./N2;
S_p=U1U2.*DU1DY+U3U2.*DU3DY;

% Compute turbulent length scales
l_O=sqrt(epsilon./N2.^(3/2));
l_b=sqrt((U1me.^2+U3me.^2)./N2);
eta=(Re^3*epsilon).^(-0.25);