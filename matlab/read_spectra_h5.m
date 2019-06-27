% Read vertical wavenumber spectra from spectra.h5

rundir='../example_run/';

read_input; Re=1/NU;

fname=[rundir 'spectra.h5'];

nk=h5readatt(fname,'/U1','Samples');

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
    
    dname=['/U1/' num];
    U1(i,:)=h5read(fname,dname);
    dname=['/U2/' num];
    U2(i,:)=h5read(fname,dname);
    dname=['/U3/' num];
    U3(i,:)=h5read(fname,dname);
    dname=['/TH1/' num];
    TH1(i,:)=h5read(fname,dname);
end

NKY=(size(U1,2)-1)/2;
E1=zeros(nk,NKY+1);     E2=zeros(nk,NKY+1);     E3=zeros(nk,NKY+1);
ETH=zeros(nk,NKY+1);

E1(:,1)=U1(:,1);      E2(:,1)=U2(:,1);
E3(:,1)=U3(:,1);      ETH(:,1)=TH1(:,1);
for i=1:NKY
    E1(:,i+1)=U1(:,i+1)+U1(:,end+1-i);
    E2(:,i+1)=U2(:,i+1)+U2(:,end+1-i);
    E3(:,i+1)=U3(:,i+1)+U3(:,end+1-i);
    ETH(:,i+1)=TH1(:,i+1)+TH1(:,end+1-i);
end
Et=E1+E2+E3;

KY=0:NKY;