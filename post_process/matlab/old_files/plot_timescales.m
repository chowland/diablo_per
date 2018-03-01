sigma_nd=0.5555/2;
sigma_d=2.5e-3;
%l_b=0.111/sqrt(300)*sigma_d/sigma_nd;
kappa_c=1e-6/300;

%v_s=[0 5 20 60 100 200 500]*10^-6;
%v_s=[0:1:600]*10^-6;
v_s=100e-6;
l_b0=0.111/sqrt(300)*sigma_d/sigma_nd;
l_b=[l_b0*5 l_b0*3 l_b0*2 l_b0];
epsilon=([1.6e-9 1.2e-8 6.25e-8 1e-6]);
l_s=1e-3;

tau_m=(4*l_b).^2/kappa_c;
tau_c=4*l_b./v_s;
tau_u=50*ones(length(epsilon),1);


%plot(v_s,tau_m,'r.-');
%hold on
%plot(v_s,tau_c,'b.-');
%plot(v_s,tau_u,'g-');

plot(epsilon,tau_m,'r.-');
hold on
plot(epsilon,tau_c,'b.-');
plot(epsilon,tau_u,'g.-');


figure
% Plot the scaling for the uptake with epsilon
v_s=100e-6;
uptake=[7.5008e-17 1.3028e-16 1.2812e-16 1.1652e-16 1.1004e-16 1.9271e-16 1.6741e-17 1.1652e-16];
%l_b=[l_b0 2*l_b0 3*l_b0 5*l_b0 51.125e-6 5*l_b0];
l_b=[l_b0 2*l_b0 3*l_b0 5*l_b0 60e-6 5*l_b0 5*l_b0 5*l_b0];
tau_c=4*l_b./v_s;
tau_m=(4*l_b).^2/kappa_c;
tau_u=[100 100 100 100 100 1000 10 100];
R1=tau_m./tau_c;
R2=tau_u./tau_c;
c=0.008207;
test=1./(R1.^-1+R2.^-1);
%test=1./(R1.^-1+R2.^-1)*0.175.*(l_s./(4*l_b));
%test=1./(R1.^-1+R2.^-1)*0.175/sqrt(300);
total=1.7e-15;
plot(test(1:5),uptake(1:5)/total,'ro')
hold on
plot(test(6:8),uptake(6:8)/total,'bs')

% Plot the scaling for the uptake with v_s
uptake=[0 4.9787e-18 1.6215e-17 4.7805e-17 7.5008e-17 1.3584e-16 2.4335e-16 3.1397e-16];
clear v_s;
v_s=[0 5 20 60 100 200 500 500]*1e-6;
tau_c=4*l_b0./v_s;
tau_m=(4*l_b0)^2/kappa_c;
tau_u=[100 100 100 100 100 100 100 15];
R1=tau_m./tau_c;
R2=tau_u./tau_c;
test=1./(R1.^-1+R2.^-1);
%test=1./(R1.^-1+R2.^-1)*0.175/sqrt(300);
%test=1./(R1.^-1+R2.^-1)*0.175.*l_b/l_s;
%test=(l_s/l_b0)./tau_c+2*l_b0*(tau_c./tau_m.^2)/l_s
plot(test(1:7),uptake(1:7)/total,'g.')
plot(test(8),uptake(8)/total,'bo');


figure
uptake=[7.5008e-17 1.3028e-16 1.2812e-16 1.1652e-16 1.1004e-16 4.9787e-18 1.6215e-17 4.7805e-17 7.5008e-17 1.3584e-16 2.4335e-16 1.927e-16];
%uptake=[7.5008e-17 1.3028e-16 1.2812e-16 1.1652e-16 1.9271e-16 4.9787e-18 1.6215e-17 4.7805e-17 1.3584e-16];
nonmotile_uptake=[8.1321e-16 7.7809e-16 7.809e-16 7.9015e-16 8.0309e-16 8.55049e-16 8.4214e-16 8.2718e-16 7.8432e-16 7.287e-16 7.4898e-16 0];
motile_uptake=[8.8822e-16 9.0837e-16 9.0902e-16 9.0668e-16 9.1313e-16 8.60024e-16 8.5836e-16 8.7498e-16 9.2016e-16 9.7206e-16 9.417e-16 3.1397e-16];
v_s=[100 100 100 100 100 5 20 60 200 500 100 500]*1e-6;
l_b=[l_b0 2*l_b0 3*l_b0 5*l_b0 60e-6 l_b0 l_b0 l_b0 l_b0 l_b0 5*l_b0 l_b0];
%tau_u=[100 100 100 100 100 100 100 100 100];
tau_u=[100 100 100 100 100 100 100 100 100 100 1000 15];
tau_c=4*l_b./v_s;
tau_m=(4*l_b).^2/kappa_c;
R1=tau_m./tau_c;
R2=tau_u./tau_c;
%test=min(R1,R2);
test=1./(R1.^-1+R2.^-1);
%test=1./(7.05e16./R1+5.74e16./R2);
%test=(motile_uptake+nonmotile_uptake)./(1-2./R2-2./R1);
%test=R2+R2;
uptake_max=1.7e-15;
y=(motile_uptake-nonmotile_uptake)/total;
c_test=(y./R1+(y.^2./R1.^2+4*y./R2).^0.5)/2;
c_test_bar=mean(c_test(6:8));
test=1./((R1.*c_test_bar.^2).^-1+(R2.*c_test_bar).^-1);
plot(test,(motile_uptake-nonmotile_uptake)/total,'r.');




