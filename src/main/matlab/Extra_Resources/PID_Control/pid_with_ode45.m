clear, close all
%
% MATLAB script file pid_with_ode45.m
%
% Written by M. Foley, January 28/21
%
global TIME_PLOT Q_PLOT FIN_PLOT FOUT_PLOT TIN_PLOT VSP_PLOT TSP_PLOT
global FSS QSS VSS TSS TINSS RHO CP TREF KC1 KC2 TI2
global TIN_STEPTIME TIN_STEP TSP_STEPTIME1 TSP_STEP1 TSP_STEPTIME2 TSP_STEP2
global FIN_STEPTIME FIN_STEP Q_MIN Q_MAX FOUT_MIN FOUT_MAX
%
TIME_PLOT = []; Q_PLOT = []; FIN_PLOT = []; FOUT_PLOT = [];
TIN_PLOT = []; VSP_PLOT = []; TSP_PLOT = [];
%
tspan = [0 45]'; % <h>
TINSS = 30; % <C>
FSS = 300; % <m3/h>
vmax = 200; % <m3>
TSS = 70; % <C>
RHO = 1000; % <kg/m3>
CP = 4.186; % <kJ/kg K>
TIN_STEPTIME = 1; % <h>
TIN_STEP = 0; % <C>
TSP_STEPTIME1 = 1; % <h>
TSP_STEP1 = 5; % <C>
TSP_STEPTIME2 = 12; % <h>
TSP_STEP2 = -10; % <C>
FIN_STEPTIME = 25; % <h>
FIN_STEP = -25; % <m3/h>
%
% Level and temperature controller tuning
%
KC1 = -2; % <(m3/h)/m3> Direct-acting, so KC1 < 0
KC2 = 1e6; % <(kJ/h)/C> Reverse-acting, so KC2 > 0
TI2 = 0.5; % <h>
%
%--------------------------------------------------------------------------
%
VSS = vmax/2; % <m3>
TREF = 0; % <C>
%
% Steady-state energy balance
%
QSS = FSS*RHO*CP*(TSS-TINSS);
FOUT_MIN = 0; FOUT_MAX = 2*FSS; % <m3/h>
Q_MIN = 0; Q_MAX = 2*QSS; % <kJ/h>
%
x0 = [VSS TSS 0]';
%
[tout,x] = ode45('pid_odefun',tspan,x0);
%
v = x(:,1); t = x(:,2); 
dummy = [TIME_PLOT Q_PLOT FIN_PLOT FOUT_PLOT TIN_PLOT VSP_PLOT TSP_PLOT];
dummy2 = sortrows(dummy,1);
time_plot = dummy2(:,1); q_plot = dummy2(:,2); 
fin_plot = dummy2(:,3); fout_plot = dummy2(:,4); 
tin_plot = dummy2(:,5); vsp_plot = dummy2(:,6); tsp_plot = dummy2(:,7); 
%
figure(1), clf
subplot(211), plot(tout,v), hold on, plot(time_plot,vsp_plot,'r'), hold off
axis([tout(1) tout(end) 70 130])
grid on, legend('V(t)','Vsp(t)')
title('Proportional Control of Liquid Volume in a Stirred-Tank Heater')
ylabel('Liquid Volume and Setpoint <m3>')
%
subplot(212), plot(time_plot,fin_plot,time_plot,fout_plot,'r'), grid on
axis([time_plot(1) time_plot(end) 270 330])
legend('Fin(t)','Fout(t)')
ylabel('Inlet & Outlet Flowrates <m3/h>'), xlabel('Time <h>')
%
figure(2), clf
subplot(311), plot(tout,t), hold on, plot(time_plot,tsp_plot,'r'), hold off
axis([tout(1) tout(end) 55 80])
grid on, legend('T(t)','Tsp(t)')
title('PI Control of Temperature in a Stirred-Tank Heater')
ylabel('Outlet Temperature <C>')
%
subplot(312), plot(time_plot,q_plot), grid on
axis([tout(1) tout(end) Q_MIN Q_MAX])
ylabel('Heat Input Rate <kJ/h>')
%
subplot(313), plot(time_plot,fin_plot,time_plot,fout_plot,'r'), grid on
axis([time_plot(1) time_plot(end) 270 330])
legend('Fin(t)','Fout(t)')
ylabel('Flowrates <m3/h>'), xlabel('Time <h>')
%
figure(1)