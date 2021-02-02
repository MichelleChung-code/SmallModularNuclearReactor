clear, close all
%
% MATLAB script file pid_with_simulink.m
%
% Written by M. Foley, January 28/21
%
global RHO CP TREF
%
tspan = [0 45]'; % <h>
tinss = 30; % <C>
fss = 300; % <m3/h>
vmax = 200; % <m3>
tss = 70; % <C>
RHO = 1000; % <kg/m3>
CP = 4.186; % <kJ/kg K>
tin_steptime = 1; % <h>
tin_step = 0; % <C>
tsp_steptime1 = 1; % <h>
tsp_step1 = 5; % <C>
tsp_steptime2 = 12; % <h>
tsp_step2 = -10; % <C>
fin_steptime = 25; % <h>
fin_step = -25; % <m3/h>
%
% Level and temperature controller tuning
%
kc1 = -2; % <(m3/h)/m3> Direct-acting, so KC1 < 0
kc2 = 1e6; % <(kJ/h)/C> Reverse-acting, so KC2 > 0
ti2 = 0.5; % <h>
%
%--------------------------------------------------------------------------
%
vss = vmax/2; % <m3>
TREF = 0; % <C>
%
% Steady-state energy balance
%
qss = fss*RHO*CP*(tss-tinss);
fout_min = 0; fout_max = 2*fss; % <m3/h>
q_min = 0; q_max = 2*qss; % <kJ/h>
%
x0 = [vss tss]';
save x0.mat x0
%
sim('pid_with_sfun');
%
figure(1), clf
subplot(211), plot(tout,v,tout,vsp_plot,'r'), hold off
axis([tout(1) tout(end) 70 130])
grid on, legend('V(t)','Vsp(t)')
title('Proportional Control of Liquid Volume in a Stirred-Tank Heater')
ylabel('Liquid Volume and Setpoint <m3>')
%
subplot(212), plot(tout,fin_plot,tout,fout_plot,'r'), grid on
axis([tout(1) tout(end) 270 330])
legend('Fin(t)','Fout(t)')
ylabel('Inlet & Outlet Flowrates <m3/h>'), xlabel('Time <h>')
%
figure(2), clf
subplot(311), plot(tout,t,tout,tsp_plot,'r'), hold off
axis([tout(1) tout(end) 55 80])
grid on, legend('T(t)','Tsp(t)')
title('PI Control of Temperature in a Stirred-Tank Heater')
ylabel('Outlet Temperature <C>')
%
subplot(312), plot(tout,q_plot), grid on
axis([tout(1) tout(end) q_min q_max])
ylabel('Heat Input Rate <kJ/h>')
%
subplot(313), plot(tout,fin_plot,tout,fout_plot,'r'), grid on
axis([tout(1) tout(end) 270 330])
legend('Fin(t)','Fout(t)')
ylabel('Flowrates <m3/h>'), xlabel('Time <h>')
%
figure(1)