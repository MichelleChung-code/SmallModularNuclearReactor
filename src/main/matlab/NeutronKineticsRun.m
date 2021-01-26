clc, clear, close all 
warning('off','all')

% Coupling Coefficients supporting the nodal approach  
coupling_coeffs_matrix = [3.5 7.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          1.8 7.1 5.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 2.2 6.9 3.9 0.0 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 3.1 7.0 3.3 0.0 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 3.7 7.0 2.9 0.0 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 4.2 7.0 2.8 0.0 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 4.5 7.0 2.6 0.0 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 4.7 7.0 2.4 0.0;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 7.0 2.2;
                          0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.6 3.5]*10^-3;

% Number of Nodes
N = 10; % cannot do 1 node yet because in the code we sometimes index with obj.N-1 which an index of 0 is invalid for MATLAB
% coupling_coeffs_matrix = compute_coupling_coefficients(N);

disp("Starting to Solve Equations");

tspan = [0 10000]';

from_csv = readtable(strcat("Initial Values", "_", string(N), "_Nodes.csv"));
csv_array = table2array(from_csv(:,4));
x0 = csv_array;


% Step change in control rod reactivity
natural_implied_control_rod_reactivity = 0.03419; 
[x0_row_num, x0_col_num] = size(x0);

control_rod_reactivity_step_size = natural_implied_control_rod_reactivity * 0.05; % 5% of natural reactivity
control_rod_reactivity_step_time = 2500; % time in seconds

% if no step response desired, just overwrite with 0, i.e. uncomment the
% line below

% control_rod_reactivity_step_size = 0; 

neutron_kinetics = NeutronKinetics(coupling_coeffs_matrix, N, control_rod_reactivity_step_size, control_rod_reactivity_step_time);
[tout, x] = neutron_kinetics.solve_neutron_kinetics(tspan, x0);
disp("Solving Completed");

% NeutronKineticsPlotting.m
run NeutronKineticsPlotting
