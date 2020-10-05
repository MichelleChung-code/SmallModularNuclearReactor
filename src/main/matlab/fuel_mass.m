function [outputArg1] = fuel_mass(d_fuel_kernel,ls_coating_layer_t, ls_coating_layer_rho)
% Calculate the fuel mass

% CONSTANTS
N_fuel_elements = 479358; %number of fuel elements in reactor core @ equilibrium state
N_coated_particles_graphite_matrix = 12000; % number of coated paricles in graphite matrix
graphite_layer_t = 5e-3; % m
graphite_matrix_d = 5e-2;  %m
graphite_rho = 1.73*100^3; % in g/(m^3)

% Calculate volumes
ls_volumes = [];
running_vol = (1/6) * pi * d_fuel_kernel^3;
running_d = d_fuel_kernel;
for i=1:length(ls_coating_layer_t)
   if i==1
      ls_volumes(i) = running_vol; 
   else
      running_d = running_d + (2*ls_coating_layer_t(i));
      new_vol = (1/6) * pi * running_d^3;
      ls_volumes(i) = new_vol - running_vol; % volume of the layer
      running_vol = new_vol;
      
   end
end

% Calculate masses
ls_masses_coated_particle = [];
for i=1:length(ls_volumes)
    ls_masses_coated_particle(i) = ls_volumes(i) * ls_coating_layer_rho(i);
end

mass_coated_particle = sum(ls_masses_coated_particle); 
vol_coated_particle = sum(ls_volumes);


% N_coated_particles_graphite_matrix in the graphite matrix 
vol_graphite_matrix = (1/6) * pi * graphite_matrix_d^3;
vol_fuel_element = (1/6) * pi * (graphite_matrix_d + 2*graphite_layer_t)^3;
vol_graphite_layer = vol_fuel_element - vol_graphite_matrix;

% In the graphite matrix 
vol_all_coated_particles = vol_coated_particle * N_coated_particles_graphite_matrix;
mass_all_coated_particles = mass_coated_particle * N_coated_particles_graphite_matrix;

% remaining volume for graphite in the matrix
vol_graphite_in_matrix = vol_graphite_matrix - vol_all_coated_particles;
mass_graphite = (vol_graphite_in_matrix + vol_graphite_layer) * graphite_rho;

disp(vol_graphite_in_matrix*12000)
outputArg1 = N_fuel_elements*(mass_graphite + mass_all_coated_particles)/1000;

% currently returns 8.9864e+04
% when:
% d_fuel_kernel = 2*0.25e-3; % in m
% ls_coating_layer_t = [0, 0.09e-3, 0.04e-3, 0.035e-3, 0.04e-3]; % in m
% ls_coating_layer_rho = [10.4, 1.1, 1.9, 3.18, 1.9];
% ls_coating_layer_rho = ls_coating_layer_rho*100^3; % in g/(m^3)
% 
% 
% fuelMass = fuel_mass(d_fuel_kernel, ls_coating_layer_t, ls_coating_layer_t) % in kg

end

