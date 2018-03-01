% Produce energy and dissipation plots mirroring those in Furue(2003)

% To be run after read_stats_h5.m
U=0.1;          % Velocity scale
L=20;           % Length scale
T=L/U;          % Time scale
V=(2*pi)^3;     % Nondimensional volume

% Formatting preamble
set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultTextboxshapeInterpreter','latex');

figure(1);
plot(T*tii,U^2*0.5*(urms+vrms+wrms+thvar/V));
hold on
plot(T*tii,U^2*0.5*(urms+vrms+wrms),'--');
plot(T*tii,U^2*0.5*thvar/V,'-.');
xlabel('$t [s]$')
ylabel('$E \ [m^2s^{-2}]$')
title('Domain-average energy')
legend('Total','Kinetic','Potential')
xlim([0 T*tii(end)])
ylim([0 3e-5])

figure(2);
plot(T*tii,U^2/T*(chi+epsilon));
hold on
plot(T*tii,U^2/T*epsilon,'--');
plot(T*tii,U^2/T*chi,'-.');
xlabel('$t [s]$')
ylabel('$\epsilon \ [m^2s^{-3}]$')
title('Dissipation rates')
legend('Total','Kinetic','Potential')
xlim([0 T*tii(end)])
ylim([0 1.5e-9])