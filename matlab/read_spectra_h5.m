% Read vertical wavenumber spectra from spectra.h5

rundir='/local/scratch/public/cjh225/test_force/Maffioli/';

read_input; Re=1/NU;

fname=[rundir 'spectra.h5'];

nk=h5readatt(fname,'/U1L','Samples');

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
    
    dname=['/U1L/' num];
    U1L(i,:)=h5read(fname,dname);
    dname=['/U1S/' num];
    U1S(i,:)=h5read(fname,dname);
    dname=['/U2L/' num];
    U2L(i,:)=h5read(fname,dname);
    dname=['/U2S/' num];
    U2S(i,:)=h5read(fname,dname);
    dname=['/U3L/' num];
    U3L(i,:)=h5read(fname,dname);
    dname=['/U3S/' num];
    U3S(i,:)=h5read(fname,dname);
    dname=['/TH1L/' num];
    TH1L(i,:)=h5read(fname,dname);
    dname=['/TH1S/' num];
    TH1S(i,:)=h5read(fname,dname);
end

NKY=(size(U1L,2)-1)/2;
E1=zeros(nk,NKY+1);     E2=zeros(nk,NKY+1);     E3=zeros(nk,NKY+1);
ETH=zeros(nk,NKY+1);

% E1(:,1)=U1(:,1);      E2(:,1)=U2(:,1);
% E3(:,1)=U3(:,1);      ETH(:,1)=TH1(:,1);
E1(:,1)=U1L(:,1)+U1S(:,1);      E2(:,1)=U2L(:,1)+U2S(:,1);
E3(:,1)=U3L(:,1)+U3S(:,1);      ETH(:,1)=TH1L(:,1)+TH1S(:,1);
for i=1:NKY
%     E1(:,i+1)=U1(:,i+1)+U1(:,end+1-i);
%     E2(:,i+1)=U2(:,i+1)+U2(:,end+1-i);
%     E3(:,i+1)=U3(:,i+1)+U3(:,end+1-i);
%     ETH(:,i+1)=TH1(:,i+1)+TH1(:,end+1-i);
    E1(:,i+1)=U1L(:,i+1)+U1L(:,end+1-i)+U1S(:,i+1)+U1S(:,end+1-i);
    E2(:,i+1)=U2L(:,i+1)+U2L(:,end+1-i)+U2S(:,i+1)+U2S(:,end+1-i);
    E3(:,i+1)=U3L(:,i+1)+U3L(:,end+1-i)+U3S(:,i+1)+U3S(:,end+1-i);
    ETH(:,i+1)=TH1L(:,i+1)+TH1L(:,end+1-i)+TH1S(:,i+1)+TH1S(:,end+1-i);
end
Et=E1+E2+E3;

E_large=zeros(nk,NKY+1);    E_small=zeros(nk,NKY+1);
E_large(:,1)=U1L(:,1)+U2L(:,1)+U3L(:,1);
E_small(:,1)=U1S(:,1)+U2S(:,1)+U3S(:,1);
for i=1:NKY
    E_large(:,i+1)=U1L(:,i+1)+U1L(:,end+1-i) ...
        +U2L(:,i+1)+U2L(:,end+1-i) + U3L(:,i+1)+U3L(:,end+1-i);
    E_small(:,i+1)=U1S(:,i+1)+U1S(:,end+1-i) ...
        +U2S(:,i+1)+U2S(:,end+1-i) + U3S(:,i+1)+U3S(:,end+1-i);
end

Etot=E_large+E_small;
KY=0:NKY;