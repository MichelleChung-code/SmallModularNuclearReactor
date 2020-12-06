clear, close all
clc
base_path = fileparts(pwd)
unit_conversion_path = strcat(base_path, '\PreliminaryReactor\UnitConversionFactors.m')
run(unit_conversion_path);


% All units are in english units
% input data is for explanatory purposes - not accurate for the specific
% model
beta = 0.0067;
beta1 = 0.000256;
beta2 = 0.0014;
beta3 = 0.0013;
beta4 = 0.0027;
beta5 = 0.00086;
beta6 = 0.00017;

mean_prompt_neutron_gen_t = 0.07669; 
lambda1 = 0.0124;
lambda2 = 0.0305;
lambda3 = 0.111;
lambda4 = 0.301;
lambda5 = 1.14;
lambda6 = 3.01;

alpha_c = (CelsiusToFahrenheit((-4.36E-5)^-1))^-1; % 1/F 
alpha_f = (CelsiusToFahrenheit((-9.40E-6)^-1))^-1;  % 1/F
theta_1o = CelsiusToFahrenheit(250); % F
theta_2o = CelsiusToFahrenheit(700); % F
Tfo = 1382;% this is 750 decrees celius 

A_fc = 4.7501e+03*m2_over_ft2; % ft^2
A_fc1 = A_fc/2; % ft^2
A_fc2 = A_fc/2; % ft^2

%http://www.endmemo.com/sconvert/j_kgcbtu_lbf.php#:~:text=F)%E2%86%94J%2F(g.C,Btu%2F(lb.
Cp_c = 1.240327; % 5.193e3 J/kgK converted to Btu/lbm F
Cp_f = 0.387277; % 1621.45 J/kgK converted to Btu/lbm F 

f = 0.97; %fuel efficiency factor
mc = 131.8847497*kg_over_lbm; % in lbm
mc1 = mc/2; 
mc2 = mc/2; 
mf =  8.5367e+04*kg_over_lbm; % in lbm
Po = 250*MWatt_over_Btu_s; % in Btu/s

% http://www.endmemo.com/convert/heat%20transfer%20coefficient.php
Ufc = 0.44583; % 9.1075e+03 W/m^2K converted to 0.44583 Btu/(s ft^2 F)

Wc = 145*kg_over_lbm; % lbm/s
T_inlet = CelsiusToFahrenheit(538); % F 

sim_path = strcat(base_path, '\PreliminaryReactor\SMR_simplified_core.slx');
disp(sim_path)
disp('Starting to run simplified reactor core')
sim(sim_path);
disp('Model has run')

% Some Plotting
figure(1)
nexttile, plot(ans.tout, ans.Tf), grid on
title('Fuel Temperature (Tf)')
xlabel('time'), ylabel('Temperature (F)')

nexttile, plot(ans.tout, ans.rho), grid on
title('Reactivity (\rho)')
xlabel('time'), ylabel('\rho')
