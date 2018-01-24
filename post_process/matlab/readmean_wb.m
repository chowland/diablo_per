% Reads in statistics outputted by diablo 
turbo=dlmread('../../control_1/mean_wb.txt');

% Enter the viscosity
NU=0.00133;

% Enter the number of scalars
N_TH=0;

% Enter the richardson number for each scalar
RI(1)=0;
RI(2)=0;
RI(3)=0;

% Set the starting time in code units for start of averaging
tstart=0;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/2)

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  E_0(k)=turbo(row,1);
  E_W(k)=turbo(row,2);
  E_S(k)=turbo(row,3);
  EPS_L(k)=turbo(row,4);
  EPS_S(k)=turbo(row,5);
  EPS_S_H(k)=turbo(row,6);
  EPS_S_V(k)=turbo(row,7);
  row=row+1;
end

