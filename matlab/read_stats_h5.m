% Read bulk quantities from stats.h5

rundir='../example_run/';

read_input; Re=1/NU;

fname=[rundir 'stats.h5'];

nk=h5readatt(fname,'/U1rms','Samples');

tii=zeros(1,nk);
U1rms=zeros(1,nk);  U2rms=zeros(1,nk);  U3rms=zeros(1,nk);
THrms=zeros(1,nk);  THflux_av=zeros(1,nk);
epsilon_av=zeros(1,nk);    chi_av=zeros(1,nk);

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

    % Read bulk quantities
    dname=['/U1rms/' num];
    tii(i)=h5readatt(fname,dname,'Time');
    U1rms(i)=h5read(fname,dname);
    dname=['/U2rms/' num];
    U2rms(i)=h5read(fname,dname);
    dname=['/U3rms/' num];
    U3rms(i)=h5read(fname,dname);
    dname=['/THrms/' num];
    THrms(i)=h5read(fname,dname);
    dname=['/THflux/' num];
    THflux_av(i)=h5read(fname,dname);
    dname=['/epsilon/' num];
    epsilon_av(i)=h5read(fname,dname);
    dname=['/chi/' num];
    chi_av(i)=h5read(fname,dname);
end

% Compute further quantities
eta_av=(Re^3*epsilon_av).^(-1/4);
EK=0.5*(U1rms.^2+U2rms.^2+U3rms.^2);
EP=0.5*Ri_t*THrms.^2;
E=EK+EP;
Re_lambda=EK.*sqrt(15*Re./epsilon_av);
Re_b_av=epsilon_av*Re/Ri_t;