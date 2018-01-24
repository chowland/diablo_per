NX=64;
NY_S=31;
LX=2*pi;
LY=2*pi;
% NX=130;
NY=64;
% Get x-vector and y-vector with uniform spacing
for i=1:NX
  xvec(i)=LX*i/NX;
end
for j=1:NY-1
  yvec(j)=LY*j/NY;
end
yvec(NY)=yvec(NY-1)+yvec(NY-1)-yvec(NY-2);

% for i=1:NX_S+1
%   xvec1(i)=xvec(i);
%   xvec2(i)=xvec(i);  
%   xvec3(i)=xvec(i+65);
%   xvec4(i)=xvec(i+65);
% end  
for j=1:NY_S+1
  yvec1(j)=yvec(j);
  yvec2(j)=yvec(j+NY_S+1);
end

mat_u=[mat_u1 mat_u3]; mat_v=[mat_u2 mat_u4];
c = 0.8/sqrt(2);

for k=1:nk
% surf(xvec1,yvec1,zeros(NX_S+1,NY_S+1)',mat_u1(1:NX_S+1,1:NY_S+1,k)','EdgeColor','none'),view(0,90)
A1 = mat_u(:,:,k)'; A2 = mat_v(:,:,k)';
subplot(1,2,1);
pcolor(xvec,yvec,A1); shading interp
caxis([-c c]); colorbar;
title('$z=0$'); xlabel('$x$'); ylabel('$y$');
subplot(1,2,2);
pcolor(xvec,yvec,A2); shading interp
caxis([-c c]); colorbar;
title('$z=\pi$'); xlabel('$x$'); ylabel('$y$');
a = annotation('textbox',[0.48 0.8 0.1 0.1],'String',['$t=' num2str(tii(k)) '$'],'FitBoxToText','on');
% hold on
% surf(xvec2,yvec2,zeros(NX_S+1,NY_S+1)',mat_u2(1:NX_S+1,1:NY_S+1,k)','EdgeColor','none'),view(0,90)
% surf(xvec3,yvec3,zeros(NX_S+1,NY_S+1)',mat_u3(1:NX_S+1,1:NY_S+1,k)','EdgeColor','none'),view(0,90)
% surf(xvec4,yvec4,zeros(NX_S+1,NY_S+1)',mat_u4(1:NX_S+1,1:NY_S+1,k)','EdgeColor','none'),view(0,90)
M(k)=getframe(gcf);
clf
end





