NX=128;
NY=65;
N_TH=0;
LX=1;
LY=1;

% Get x-vector and y-vector with uniform spacing
for i=1:NX
  xvec(i)=LX*i/NX;
end
for j=1:NY-1
  yvec(j)=LY*j/NY;
end
yvec(NY)=yvec(NY-1)+yvec(NY-1)-yvec(NY-2);

null(1:NX,1:NY)=0;

colormat=get(gca,'ColorOrder');

moviefile=dlmread('movie_xz_3.txt');

% Get the number of timesteps
nk=length(moviefile(:,1))/(NX*NY);
%mean_size=size(ume);
%nk=mean_size(1,2);


count=0;

% zero matrix mat
mat(1:NX,1:NY,1:nk)=0;


for k=1:nk

for j=1:NY
for i=1:NX
count=count+1;
  xvec(i)=moviefile(count,1);
  yvec(j)=moviefile(count,2);
  mat(i,j,k)=moviefile(count,3);
end
end 

end

surf(xvec,yvec,zeros(NX,NY)',mat(:,:,nk)','EdgeColor','none'),view(0,90)




