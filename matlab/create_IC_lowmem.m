% Create initial condition when variables too large to all be stored in
% memory at once

% Dimensionless parameters
NU=1e-4;    Re=1/NU;    Ri_t=1;     Pr=1;
% Dimensional buoyancy frequency and viscosity
N=1e-2;     nu=1e-6;
% Dimensional length and velocity scales
L=sqrt(nu*Re*Ri_t^0.5/N);   U=sqrt(nu*Re*N/Ri_t^0.5);
% Rescale dimensional parameters for spectrum
b=1300/L;   N0=5.2e-3*L/U;      f=7.3e-5*L/U;
% Cutoff wavenumbers for shear modes & wave spectrum
KY_C=7;     K_C=7;
% Calculate constant pre-factor (for wave components)
E=6.3e-5;   Sigma=0.468;        A=Ri_t*b*E*f*sqrt(Ri_t-f^2)/(2*pi*Sigma);

N=1024;       N_CORES=256;      NK=floor(N/3);      K=[0:N/2 -(N/2-1):-1];  
NX=N; NY=N; NZ=N; NKX=NK; NKY=NK; NKZ=NK;
KX=K; KY=K; KZ=K; kappa=zeros(NX,NZ);
for i=1:NX; for k=1:NZ
        kappa(i,k)=sqrt(KX(i)^2+KZ(k)^2);
end; end

CS1=zeros(NX,NY,NZ);

% Define unique values (truncated to wavenumbers NK)
for i=2:NKX+1
    for k=[1:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
            K2=kappa(i,k)^2+KY(j)^2;
            if K2<=K_C^2
                alpha=2*pi*rand;
                CS1(i,j,k)=sqrt(A*KY(j)^2/kappa(i,k)/ ...
                    sqrt(K2)/(Ri_t*kappa(i,k)^2+f^2*KY(j)^2) ...
                    /(KY(j)^2+pi^2*Ri_t*9/(b^2*N0^2)))*exp(1i*alpha);
            end
        end
    end
end
for k=2:NKZ+1
    for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
        K2=kappa(1,k)^2+KY(j)^2;
        if K2<=K_C^2
            alpha=2*pi*rand;
            CS1(1,j,k)=sqrt(A*KY(j)^2/kappa(1,k)/ ...
                sqrt(K2)/(Ri_t*kappa(1,k)^2+f^2*KY(j)^2) ...
                /(KY(j)^2+pi^2*Ri_t*9/(b^2*N0^2)))*exp(1i*alpha);
        end
    end
end
for j=2:NKY+1      % SHEAR FLOW COMPONENT DEFINED HERE
    if KY(j)^2<=KY_C^2      % Set so that mean Ri_g = Ri_tau
    CS1(1,j,1)=1./(sqrt(2*KY_C)*KY(j));
    beta(j)=rand;
    end
end
CS1(1,1,1)=0;
disp('CS1 calculated')
ES=sum(abs(CS1(1,:,1)).^2);
EK=sum(abs(CS1(:)).^2)-ES;
A0=sqrt(ES/18/EK);

%% U1
CU1=zeros(NX,NY,NZ);
for i=2:NKX+1
    for k=[1:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
            K2=kappa(i,k)^2+KY(j)^2;
            if K2<=K_C^2
                CU1(i,j,k)=CS1(i,j,k)*KY(j)*KX(i)/sqrt(K2) ...
                    /kappa(i,k);
            end
        end
    end
end
for k=2:NKZ+1
    for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
        K2=kappa(1,k)^2+KY(j)^2;
        if K2<=K_C^2
            CU1(1,j,k)=CS1(1,j,k)*KY(j)*KX(1)/sqrt(K2) ...
                /kappa(1,k);
        end
    end
end
for j=2:NKY+1      % SHEAR FLOW COMPONENT DEFINED HERE
    if KY(j)^2<=KY_C^2      % Set so that mean Ri_g = Ri_tau
    alpha=2*pi*rand;
    CU1(1,j,1)=sqrt(beta(j))*CS1(1,j,1)*exp(1i*alpha);
    end
end
CU1=A0*CU1;
CU1(1,:,1)=CU1(1,:,1)/A0;
for i=NX-NKX+2:NX
    for k=[2:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY]
            CU1(i,j,k)=conj(CU1(NX-i+2,NY-j+2,NZ-k+2));
            end
    end
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU1(i,j,1)=conj(CU1(NX-i+2,NY-j+2,1));
    end
end
for k=NZ-NKZ+2:NZ
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU1(1,j,k)=conj(CU1(1,NY-j+2,NZ-k+2));
    end
end
for j=NY-NKY+2:NY
    CU1(1,j,1)=conj(CU1(1,NY-j+2,1));
end
U1=ifftn(CU1*NX*NY*NZ,'symmetric'); clear CU1; disp('Transformed U1');
h5create('start.h5','/U1',[NX NY NZ],'ChunkSize',[NX 1 NZ/sqrt(N_CORES)])
h5write('start.h5','/U1',U1)
clear U1

%% U2
CU2=zeros(NX,NY,NZ);
for i=2:NKX+1
    for k=[1:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
            K2=kappa(i,k)^2+KY(j)^2;
            if K2<=K_C^2
                CU2(i,j,k)=-CS1(i,j,k)*kappa(i,k)/sqrt(K2);
            end
        end
    end
end
for k=2:NKZ+1
    for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
        K2=kappa(1,k)^2+KY(j)^2;
        if K2<=K_C^2
            CU2(1,j,k)=-CS1(1,j,k)*kappa(1,k)/sqrt(K2);
        end
    end
end
CU2=A0*CU2;
for i=NX-NKX+2:NX
    for k=[2:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY]
            CU2(i,j,k)=conj(CU2(NX-i+2,NY-j+2,NZ-k+2));
            end
    end
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU2(i,j,1)=conj(CU2(NX-i+2,NY-j+2,1));
    end
end
for k=NZ-NKZ+2:NZ
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU2(1,j,k)=conj(CU2(1,NY-j+2,NZ-k+2));
    end
end
for j=NY-NKY+2:NY
    CU2(1,j,1)=conj(CU2(1,NY-j+2,1));
end
U2=ifftn(CU2*NX*NY*NZ,'symmetric'); clear CU2; disp('Transformed U2');
h5create('start.h5','/U2',[NX NY NZ],'ChunkSize',[NX 1 NZ/sqrt(N_CORES)])
h5write('start.h5','/U2',U2)
clear U2

%% U3
CU3=zeros(NX,NY,NZ);
for i=2:NKX+1
    for k=[1:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
            K2=kappa(i,k)^2+KY(j)^2;
            if K2<=K_C^2
                CU3(i,j,k)=CS1(i,j,k)*KY(j)*KZ(k)/sqrt(K2) ...
                    /kappa(i,k);
            end
        end
    end
end
for k=2:NKZ+1
    for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
        K2=kappa(1,k)^2+KY(j)^2;
        if K2<=K_C^2
            CU3(1,j,k)=CS1(1,j,k)*KY(j)*KZ(k)/sqrt(K2) ...
                /kappa(1,k);
        end
    end
end
for j=2:NKY+1      % SHEAR FLOW COMPONENT DEFINED HERE
    if KY(j)^2<=KY_C^2      % Set so that mean Ri_g = Ri_tau
        alpha=2*pi*rand;
        CU3(1,j,1)=sqrt(1-beta(j))*CS1(1,j,1)*exp(1i*alpha);
    end
end
CU3=A0*CU3;
CU3(1,:,1)=CU3(1,:,1)/A0;
for i=NX-NKX+2:NX
    for k=[2:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY]
            CU3(i,j,k)=conj(CU3(NX-i+2,NY-j+2,NZ-k+2));
            end
    end
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU3(i,j,1)=conj(CU3(NX-i+2,NY-j+2,1));
    end
end
for k=NZ-NKZ+2:NZ
    for j=[2:NKY+1 NY-NKY+2:NY]
        CU3(1,j,k)=conj(CU3(1,NY-j+2,NZ-k+2));
    end
end
for j=NY-NKY+2:NY
    CU3(1,j,1)=conj(CU3(1,NY-j+2,1));
end
U3=ifftn(CU3*NX*NY*NZ,'symmetric'); clear CU3; disp('Transformed U3');
h5create('start.h5','/U3',[NX NY NZ],'ChunkSize',[NX 1 NZ/sqrt(N_CORES)])
h5write('start.h5','/U3',U3)
clear U3

%% TH
CTH=zeros(NX,NY,NZ);
for i=2:NKX+1
    for k=[1:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
            K2=kappa(i,k)^2+KY(j)^2;
            if K2<=K_C^2
                CTH(i,j,k)=CS1(i,j,k)/sqrt(Ri_t)*1i;
            end
        end
    end
end
for k=2:NKZ+1
    for j=[2:NKY+1 NY-NKY+2:NY] % NB general spectrum index would start at 1
        K2=kappa(1,k)^2+KY(j)^2;
        if K2<=K_C^2
            CTH(1,j,k)=CS1(1,j,k)/sqrt(Ri_t)*1i;
        end
    end
end
CTH=A0*CTH;
for i=NX-NKX+2:NX
    for k=[2:NKZ+1 NZ-NKZ+2:NZ]
        for j=[2:NKY+1 NY-NKY+2:NY]
            CTH(i,j,k)=conj(CTH(NX-i+2,NY-j+2,NZ-k+2));
            end
    end
    for j=[2:NKY+1 NY-NKY+2:NY]
        CTH(i,j,1)=conj(CTH(NX-i+2,NY-j+2,1));
    end
end
for k=NZ-NKZ+2:NZ
    for j=[2:NKY+1 NY-NKY+2:NY]
        CTH(1,j,k)=conj(CTH(1,NY-j+2,NZ-k+2));
    end
end
for j=NY-NKY+2:NY
    CTH(1,j,1)=conj(CTH(1,NY-j+2,1));
end
TH=ifftn(CTH*NX*NY*NZ,'symmetric'); clear CTH; disp('Transformed TH');
h5create('start.h5','/TH1',[NX NY NZ],'ChunkSize',[NX 1 NZ/sqrt(N_CORES)])
h5write('start.h5','/TH1',TH)
clear TH

h5writeatt('start.h5','/','Resolution',[NX NY NZ; NX NY NZ]')
h5writeatt('start.h5','/','Date',datestr(now))
h5writeatt('start.h5','/','Info','GM spectrum created in Matlab')
h5writeatt('start.h5','/','Time',double(0))
h5writeatt('start.h5','/','Timestep',0)