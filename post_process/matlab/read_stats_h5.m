% Run directory
rundir = '/local/scratch/public/cjh225/force_run/';
% Movie file
fname = [rundir 'stats.h5'];

% Number of samples
nk = h5readatt(fname,'/urms','Samples');


% Preallocation
tii=zeros(1,nk); urms=zeros(1,nk); vrms=zeros(1,nk); wrms=zeros(1,nk);
thme=zeros(1,nk); thth=zeros(1,nk); thvar=zeros(1,nk);
epsilon=zeros(1,nk); eta=zeros(1,nk); re_lambda=zeros(1,nk);

for i=1:nk
    if i<=10
        num=['000' int2str(i-1)];
    elseif i<=100
        num=['00' int2str(i-1)];
    elseif i<=1000
        num=['0' int2str(i-1)];
    else
        num=int2str(i-1);
    end
    dname=['/urms/' num];
    tii(i)=h5readatt(fname,dname,'Time');
    urms(i)=h5read(fname,dname);
    dname=['/vrms/' num];
    vrms(i)=h5read(fname,dname);
    dname=['/wrms/' num];
    wrms(i)=h5read(fname,dname);
    dname=['/thme/' num];
    thme(i)=h5read(fname,dname);
    dname=['/thth/' num];
    thth(i)=h5read(fname,dname);
    dname=['/thvar/' num];
    thvar(i)=h5read(fname,dname);
    dname=['/epsilon/' num];
    epsilon(i)=h5read(fname,dname);
    dname=['/eta/' num];
    eta(i)=h5read(fname,dname);
    dname=['/re_lambda/' num];
    re_lambda(i)=h5read(fname,dname);
end