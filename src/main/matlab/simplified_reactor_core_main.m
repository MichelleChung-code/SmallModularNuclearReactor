clear, close all
clc
beta = 0.0067;
beta1 = 0.000256;
beta2 = 0.0014;
beta3 = 0.0013;
beta4 = 0.0027;
beta5 = 0.00086;
beta6 = 0.00017;

mean_prompt_neutron_gen_t = 7.66e-4; %calculate (1/beta sum(betai/lambdai))^-1
lambda1 = 0.0124;
lambda2 = 0.0305;
lambda3 = 0.111;
lambda4 = 0.301;
lambda5 = 1.14;
lambda6 = 3.01;

alpha_c = -5.22e-6; 
alpha_f = -2.42e-5; 
theta_1o = 482; % need to check
theta_2o = 1382; % need to check
Tfo = 1382;% need to check

A_fc = 14120.6; % need to calculate 
A_fc1 = 7060.3; % need to calculate 
A_fc2 = 7060.3; % need to calculate 
Cp_c = 1.24;
Cp_f = 0.059; % need to find 
f = 0.97; % need to find (fuel efficiency factor?) 
mc = 6898.3; % need to find 
mc1 = 3449.1; % need to find 
mc2 = 3449.1; % need to find 
mf = 198116.207; % in lbm
Po = 434100241.3585974;
Ufc = 0.0909; % need to find
Wc = 388.014;
T_inlet = 1000.4;

base_path = fileparts(pwd);
sim_path = strcat(base_path, '\simulink\SMR_simplified_core.slx');
disp(sim_path)
disp('Starting to run simplified reactor core')
sim(sim_path);
disp('Model has run')

% Some Plotting
figure(1), subplot(3,3,1)
plot(ans.tout, ans.T_outlet), grid on
title('Toutlet')
xlabel('time(s)'), ylabel('Temperature (F)')

subplot(3,3,2), plot(ans.tout, ans.Pth), grid on
title('Fractional Reactor Power (Pth)')
xlabel('time(s)'), ylabel('Pth')

subplot(3,3,3), plot(ans.tout, ans.Tf), grid on
title('Fuel Temperature (Tf)')
xlabel('time(s)'), ylabel('Temperature (F)')

subplot(3,3,4), plot(ans.tout, ans.rho), grid on
title('Reactivity (rho)')
xlabel('time(s)'), ylabel('rho')

subplot(3,3,5), plot(ans.tout, ans.theta_1), grid on
title('Theta 1')
xlabel('time(s)'), ylabel('Temperature (F)')

subplot(3,3,6), plot(ans.tout, ans.theta_2), grid on
title('Theta 2')
xlabel('time(s)'), ylabel('Temperature (F)')