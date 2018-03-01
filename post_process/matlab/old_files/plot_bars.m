
phi=5e-12;
lp=1e-3;
sigma=lp/2;
Cd=phi/((2*pi*sigma^2)^(3/2));
Bd=5e11;
uptake=0.01/Bd;


%x1=[1e-6 6.25e-8 1.2e-8 1.6e-9 0];
x1=[0 -1 -2 -3 -4];
%x1=[1e-6 6.25e-8 1.2e-8 1.6e-9 0];
y=100;
delta_y=5;
delta_x=[0.05 0.05 0.05 0.05 0.05];
%delta_x=x1/10;
z1=[7.5008e-17 1.3028e-16 1.2812e-16 1.1652e-16 1.1004e-16];
%z1=[1.48 2.55 2.46 2.18]*uptake*Cd;
%z1=[1.48 2.55 2.46 2.18 2.01]*uptake*Cd;
%z1=[1.48 2.56 2.53 2.29 2.16]*uptake*Cd;

i=1;
fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i)-delta_x(i) x1(i)-delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i)],[z1(i) z1(i) z1(i) z1(i)],[1 0 0]);  
hold on
for i=2:length(x1);
%top cap
  fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i)-delta_x(i) x1(i)-delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i)],[z1(i) z1(i) z1(i) z1(i)],[1 0 0]);  
%sides
  fill3([y-delta_y y-delta_y y-delta_y y-delta_y],[x1(i)-delta_x(i) x1(i)-delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i)],[0 z1(i) z1(i) 0],[1 0 0]);  
%  fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i)-delta_x(i) x1(i)-delta_x(i) x1(i)-delta_x(i) x1(i)-delta_x(i)],[0 0 z1(i) z1(i)],[1 0 0]);  
  fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i)+delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i)],[0 0 z1(i) z1(i)],[1 0 0]);  
% connecting line
%  fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i-1)+delta_x(i) x1(i-1)+delta_x(i) x1(i)-delta_x(i) x1(i)-delta_x(i)],[z1(i-1) z1(i-1) z1(i) z1(i)],[1 0 0]);  
  fill3([y-delta_y y+delta_y y+delta_y y-delta_y],[x1(i-1)-delta_x(i) x1(i-1)-delta_x(i) x1(i)+delta_x(i) x1(i)+delta_x(i)],[z1(i-1) z1(i-1) z1(i) z1(i)],[1 0 0]);  
end

x2=[0 5 20 60 100 200 500];
%z2=[0 0.32 1.48 4.74]*uptake*Cd;
z2=[0 4.9787e-18 1.6215e-17 4.7805e-17 7.5008e-17 1.3584e-16 2.4335e-16];
%y=0;
y=10^-6;
delta_y=0.05;
%delta_y=10^-6/10;
clear delta_x
delta_x=5;

i=1;
fill3([x2(i)-delta_x x2(i)-delta_x x2(i)+delta_x x2(i)+delta_x],[y-delta_y y+delta_y y+delta_y y-delta_y],[z2(i) z2(i) z2(i) z2(i)],[1 0 0]);  
for i=2:length(x2);
% top cap
  fill3([x2(i)-delta_x x2(i)-delta_x x2(i)+delta_x x2(i)+delta_x],[y-delta_y y+delta_y y+delta_y y-delta_y],[z2(i) z2(i) z2(i) z2(i)],[1 0 0]);  
% sides
  fill3([x2(i)-delta_x x2(i)-delta_x x2(i)-delta_x x2(i)-delta_x],[y-delta_y y+delta_y y+delta_y y-delta_y],[0 0 z2(i) z2(i)],[1 0 0]);  
%  fill3([x2(i)-delta_x x2(i)-delta_x x2(i)+delta_x x2(i)+delta_x],[y-delta_y y-delta_y y-delta_y y-delta_y],[0 z2(i) z2(i) 0],[1 0 0]);  
  fill3([x2(i)-delta_x x2(i)-delta_x x2(i)+delta_x x2(i)+delta_x],[y+delta_y y+delta_y y+delta_y y+delta_y],[0 z2(i) z2(i) 0],[1 0 0]);  
% connecting line
  fill3([x2(i-1)+delta_x x2(i-1)+delta_x x2(i)-delta_x x2(i)-delta_x],[y-delta_y y+delta_y y+delta_y y-delta_y],[z2(i-1) z2(i-1) z2(i) z2(i)],[1 0 0]);  
end

grid on

%h=line([0 500],[-delta_y -delta_y],[0 0.0098*uptake*Cd*500]);
%set(h,'LineWidth',2);
%set(h,'Color','k');


delta_vs=2;
V_s=(0:2*delta_vs:300);
%expenditure=1.138888e-11*((V_s*2*1e-6).^2)*120;
expenditure=9.2971e-11*((V_s*1e-6).^2)*120;
%V_s=(0:2*delta_vs:300);
%expenditure=9.9574e-19*(0.3*V_s);
for i=2:length(V_s)
h=fill3([V_s(i-1) V_s(i-1) V_s(i) V_s(i)],[y-delta_y y+delta_y y+delta_y y-delta_y],[expenditure(i-1) expenditure(i-1) expenditure(i) expenditure(i)],[0 0 1]);  
set(h,'EdgeColor','none');
end
%plot3(V_s,zeros(length(V_s))-delta_y,expenditure,'k-','LineWidth',2);



axis([-10 600 -4.2 0.2 0 2.5e-16]);
%axis([-10 600 10^-9 10^-6 0 2.5e-16]);
set(gca,'FontName','Times');
set(gca,'FontSize',14);
xlabel('V_s (\mu m/s)');
ylabel('\epsilon (W/kg)');

