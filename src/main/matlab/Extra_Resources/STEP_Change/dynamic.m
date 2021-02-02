clear, clear global

global DIST_PAR PLOT_REFLUX TIME_REFLUX

DIST_PAR(1) = 1.5; % relative volatility
DIST_PAR(2) = 41; % total number of stages
DIST_PAR(3) =21; % feed stage 
DIST_PAR(4) = 1;  % feed flowrate 
DIST_PAR(5) = 0.5; % feed composition, light comp
DIST_PAR(6) = 1;% feed quality (1 = sat'd liqd,0 = sat'd vapor) 
DIST_PAR(7) = 2.706; % reflux flowrate
DIST_PAR(8) = 3.206; % reboiler vapor flowrate
DIST_PAR(9) = 5; % distillate molar hold-up
DIST_PAR(10) = 5; % bottoms molar hold-up 
DIST_PAR(11) = 0.5; % stage molar hold-up
DIST_PAR(12) = 0.01*DIST_PAR(7); % size of first reflux step change
DIST_PAR(13) = 50; % time of first reflux step change 
DIST_PAR(14) = 0; % magnitude step in vapor 
DIST_PAR(15) = 0; % time of vapor step change 
DIST_PAR(16) = 0; % magnitude of feed comp change 
DIST_PAR(17) = 0; % time of feed comp change 
DIST_PAR(18) = 0; % magnitude of feed comp change 
DIST_PAR(19) = 0; % time of feed flow change 
DIST_PAR(20) = -0.01*DIST_PAR(7); % size of 2nd reflux step change 
DIST_PAR(21) = 200; % time of second reflux step change 


ns = DIST_PAR(2);
TIME_REFLUX = []; PLOT_REFLUX = [];

load x0; % Initial guess
tspan =0.0:5:400; % time length of simulation
[tup,xkup] = ode45('dist_dyn',tspan,x0);

z = [TIME_REFLUX PLOT_REFLUX]; zz = sortrows(z,1);
time_reflux = zz(:,1); plot_reflux = zz(:,2);

DIST_PAR(12) = -DIST_PAR(12);
[tdown,xkdown] = ode45('dist_dyn',tspan,x0);

save dynamic

figure(1), clf
subplot(2,1,1), plot(tup,xkup(:,1),tdown,xkdown(:,1),'r'), grid on
title('Dynamic Response of Simulated Binary Distillation Column'),
ylabel('Mole Frxn LK in Distillate (xd)')
axis([tspan(1) tspan(end) 0.94 1])
legend('+1% Reflux','-1% Reflux')

subplot(2,1,2), plot(tup,xkup(:,ns),tdown,xkdown(:,ns),'r'), grid on
ylabel('Mole Frxn LK in Bottoms (xb)')
xlabel('Time <min>'), 
axis([tspan(1) tspan(end) 0 0.06])
% 
figure(2), clf,
plot(time_reflux,plot_reflux), grid on
title('Just so you can see that we did have two step changes in the reflux flowrate')
ylabel('Reflux flow <mol/min>'), xlabel('Time <min>')
  