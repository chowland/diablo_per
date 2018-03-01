% Reads in statistics outputted by diablo

turbo=dlmread('tke.txt');

% Set the domain size
NY=256;
NYM=NY-1;

% Set the starting time in code units for start of averaging
tstart=1;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1)))/2

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  epsilon(k)=turbo(row,1); 
  eta(k)=turbo(row,2);   
  epsilon_pseudo(k)=turbo(row,3);
  epsilon_test(k)=turbo(row,4);
%  re_lambda(k)=turbo(row,3);
%  re_lambda(k)=turbo(row,4);
  row=row+1;
end

%for k=1:nk
%  re_lambda(k)=(2/3)*tke(k)*sqrt(15/(epsilon(k)*NU));
%end
