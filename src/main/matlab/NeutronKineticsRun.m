clc, clear, close all 
warning('off','all')

% Dynamic Inputs that can change  
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

disp("Starting to Solve Equations");

tspan = [0 10000]';
N = 10; % number of nodes


from_csv = readtable('Initial Values.csv');
csv_array = table2array(from_csv(:,4));
x0 = csv_array;
step_size = 0.1;

neutron_kinetics = NeutronKinetics(coupling_coeffs_matrix, N);
[tout, x] = neutron_kinetics.solve_neutron_kinetics(tspan, x0, step_size);

disp("Solving Completed");

% NeutronKineticsPlotting.m
run NeutronKineticsPlotting
