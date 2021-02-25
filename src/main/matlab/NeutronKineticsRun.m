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
N = 10; % cannot do 1 node yet because in the code we often index with obj.N-1 resulting in an invalid index of 0
% coupling_coeffs_matrix = compute_coupling_coefficients(N);

disp("Starting to Solve Equations");

tspan = [0 10000]';

from_csv = readtable(strcat("Initial Values", "_", string(N), "_Nodes.csv"));
csv_array = table2array(from_csv(:,4));
x0 = csv_array;

natural_reactivity = 0.03419; 
reactivity_profile = [-0.0545030139935414 0.329809379569424 0.614779449823997 0.76941928451865 0.90915748167749 1.05198549992605 1.20415573523475 1.36765908626817 1.54073311718303 1.56910946268381];

% apply burnup effect
reactivity_profile = fuel_burnup(80000, N, reactivity_profile);
% reactivity_profile = [0.04	0.397849462	0.655913978	0.774193548	0.860215054	0.924731183	0.967741935	0.989247312	1	0.892473118];
%reactivity_profile = [0.6	0.6	0.655913978	0.774193548	0.860215054	0.924731183	0.967741935	0.989247312	1	0.892473118];
control_rod_reactivity = 0.03;

x0 = [x0; 0]; % append the integ val for PI controller

%For step change in reactivity
reactivity_step_size = natural_reactivity * .05; % 5% of natural reactivity
reactivity_step_time = 2500; % time in seconds
control_rod_sp_step_time = 5000; % time in seconds

% if no step response desired, just overwrite with 0, i.e. uncomment the
% line below

%reactivity_step_size = 0; 

neutron_kinetics = NeutronKinetics(coupling_coeffs_matrix, N, reactivity_step_size, reactivity_step_time, natural_reactivity,reactivity_profile,control_rod_reactivity, x0, control_rod_sp_step_time);
[tout, x] = neutron_kinetics.solve_neutron_kinetics(tspan, x0);
disp("Solving Completed");

% NeutronKineticsPlotting.m
run NeutronKineticsPlotting
