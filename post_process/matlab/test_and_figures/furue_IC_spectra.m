% Furue initial condition and Froude spectra
% Case 0 or 1r
NKX=47; NKY=47; NKZ=427;
b=1300;
N=5.2E-3;
E=6.3e-5;
S=0.468;
f=7.33e-5;
C=pi/b;
DKX=2*pi/100; DKY=2*pi/100; DKZ=2*pi/128.1;
h=sqrt((DKX/2)^2+(DKY/2)^2);
KX=(0:NKX)*DKX; KY=(0:NKY)*DKY; KZ=(0:NKZ)*DKZ;
for i=1:NKX+1
    for j=1:NKY+1
        KH(i,j)=sqrt(KX(i)^2+KY(j)^2);
    end
end
for i=1:NKX+1
    for j=1:NKY+1
        for k=1:NKZ+1
            if KZ(k)==0
                Fe(i,j,k)=0;
            elseif KH(i,j)==0 % Shear Flow component
%                 Fe(i,j,k)=DKX*DKY/(pi*h^2)*2*C*b^2*N^2*E/pi/S/(KZ(k)^2+9*C^2)*atan(sqrt(N^2-f^2)*h/f/sqrt(h^2+KZ(k)^2));
                Fe(i,j,k)=0;
            else
                Fe(i,j,k)=1/(2*pi*KH(i,j))*C*b^2*N^2*E*f*sqrt(N^2-f^2)/(pi*S)*KZ(k)^2/sqrt(KH(i,j)^2+KZ(k)^2)/(N^2*KH(i,j)^2+f^2*KZ(k)^2)/(KZ(k)^2+9*C^2);
            end
            Shear(i,j,k)=KZ(k)^2*Fe(i,j,k)/(KH(i,j)^2+KZ(k)^2);
        end
    end
end

Fr=permute(sum(sum(Shear,1),2),[1 3 2]);
Fr=KZ.^2.*Fr/N^2*DKX*DKY;
loglog(KZ,Fr)
figure(2)
E2=permute(sum(sum(Fe,1),2),[1 3 2])*DKX*DKY;
loglog(KZ,E2)

%% 

NK=83;
b0=65;
RI=1;
f0=0.014;
C0=pi/b0;
DKH=1/32;
KX=(0:NK)*DKH; KY=(0:NK)*DKH; KZ=0:NK;
for i=1:NK+1
    for j=1:NK+1
        KH(i,j)=sqrt(KX(i)^2+KY(j)^2);
    end
end
for i=1:NK+1
    for j=1:NK+1
        for k=1:NK+1
            if KZ(k)==0 || KH(i,j)==0
                Fe0(i,j,k)=0;
            else
                Fe0(i,j,k)=1/(2*pi*KH(i,j))*C0*b0^2*RI*E*f0*sqrt(RI-f0^2)/(pi*S)*KZ(k)^2/sqrt(KH(i,j)^2+KZ(k)^2)/(RI*KH(i,j)^2+f0^2*KZ(k)^2)/(KZ(k)^2+9*C0^2);
            end
            Shear0(i,j,k)=KZ(k)^2*Fe0(i,j,k)/(KH(i,j)^2+KZ(k)^2);
        end
    end
end
figure(1); hold on
Fr0=permute(sum(sum(Shear0,1),2),[1 3 2])*DKH^2;
Fr0=KZ.^2.*Fr0/RI;
loglog(KZ/20,Fr0*20)
figure(2); hold on
E0=permute(sum(sum(Fe0,1),2),[1 3 2])*DKH^2;
loglog(KZ/20,E0*0.2)