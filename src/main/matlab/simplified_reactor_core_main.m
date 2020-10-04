clear, close all
clc
beta = 0.67;
beta1 = 0.0256;
beta2 = 0.14;
beta3 = 0.13;
beta4 = 0.27;
beta5 = 0.086;
beta6 = 0.017;

mean_prompt_neutron_gen_t = 7.66e-4;
lambda1 = 0.0256;
lambda2 = 0.14;
lambda3 = 0.13;
lambda4 = 0.27;
lambda5 = 0.086;
lambda6 = 0.017;

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
f = 0.97; % need to find 
mc = 6898.3; % need to find 
mc1 = 3449.1; % need to find 
mc2 = 3449.1; % need to find 
mf = 53850.5; % need to find 
Po = 434100241.3585974;
Ufc = 0.0909; % need to find
Wc = 388.014;
T_inlet = 1000.4;

base_path = fileparts(pwd);
sim_path = strcat(base_path, '\simulink\SMR_simplified_core.slx');
disp('Starting to run simplified reactor core')
sim(sim_path);
disp('Model has run')

% figure(1), plot(tout, T_outlet), grid on
% title('Toutlet')
% xlabel('time(s)'), ylabel('Temperature (F)')
% 
% figure(2), plot(tout, Pth), grid on
% title('Fractional Reactor Power (Pth)')
% xlabel('time(s)'), ylabel('Pth')
% 
% figure(3), plot(tout, Tf), grid on
% title('Fuel Temperature (Tf)')
% xlabel('time(s)'), ylabel('Temperature (F)')
% 
% figure(4), plot(tout, rho), grid on
% title('Reactivity (rho)')
% xlabel('time(s)'), ylabel('rho')