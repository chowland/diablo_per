NX=64;
NY_S=31;
NY_mpi=NY_S+1;
N_TH=0;
LX=2*pi;
LY=pi;

rundir='../../control_1/';

% Get x-vector and y-vector with uniform spacing
for i=1:NX
  xvec(i)=LX*i/NX;
end
for j=1:NY_mpi-1
  yvec(j)=LY*j/NY_mpi;
end
yvec(NY_mpi)=yvec(NY_mpi-1)+yvec(NY_mpi-1)-yvec(NY_mpi-2);

null(1:NX,1:NY_mpi)=0;

colormat=get(gca,'ColorOrder');

moviefile=dlmread([rundir 'movie_xy_4.txt']);

% Get the number of timesteps
nk=length(moviefile(:,1))/(NX*NY_mpi);
%mean_size=size(ume);
%nk=mean_size(1,2);


count=0;

% zero matrix mat
mat_u4(1:NX,1:NY_mpi,1:nk)=0;


for k=1:nk

for j=1:NY_mpi
for i=1:NX
count=count+1;
  mat_u4(i,j,k)=moviefile(count,3);
end
end 

end

for k=1:nk
% surf(xvec,yvec,zeros(NX,NY)',mat_u1(:,:,k)','EdgeColor','none'),view(0,90)
colormap(cmocean('balance'));
pcolor(xvec,yvec,mat_u4(:,:,k)'); shading interp;
colorbar; caxis(0.8/sqrt(2)*[-1 1])
M(k)=getframe(gcf);
% clf
end




