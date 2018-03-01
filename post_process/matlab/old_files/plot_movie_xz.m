NX=128;
NY=128;
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

moviefile=dlmread('movie_xz.txt');

% Get the number of timesteps
nk=length(moviefile(:,1))/(NX*NY);
%mean_size=size(ume);
%nk=mean_size(1,2);


count=0;

% zero matrix mat
mat(1:NX,1:NY,1:nk)=0;


for k=1:nk

for i=1:NX
for j=1:NY
count=count+1;
  mat(i,j,k)=moviefile(count);
end
end 

end

for k=1:nk
surf(xvec,yvec,zeros(NX,NY)',mat(:,:,k)','EdgeColor','none'),view(0,90)
M(k)=getframe(gcf);
clf
end




