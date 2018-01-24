% Reads in statistics outputted by diablo

turbo=dlmread('tke.txt');

% Set the domain size
NY=256;
NYM=NY-1;
N_TH=0;

% Set the starting time in code units for start of averaging
tstart=1;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1)))/(2+N_TH)

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  epsilon(k)=turbo(row,1); 
  eta(k)=turbo(row,2);   
  re_lambda(k)=turbo(row,3);   
%  epsilon_pseudo(k)=turbo(row,4);
  row=row+1;
  for n=1:N_TH
    thv(k,n)=turbo(row,1);
    row=row+1;
  end
end

%for k=2:nk-1
%  dtkedt(k)=(tke(k+1)-tke(k-1))/(tii(k+1)-tii(k-1));
%end

plot(tii(3:nk-1),dtkedt(3:nk-1),'k-');
hold on
plot(tii(3:nk-1),thv(3:nk-1,1),'r-');
plot(tii(3:nk-1),-epsilon(3:nk-1),'b-');
plot(tii(3:nk-1),-epsilon_pseudo(3:nk-1),'b--');
plot(tii(3:nk-1),dtkedt(3:nk-1)-squeeze(thv(3:nk-1,1))','k.-');

