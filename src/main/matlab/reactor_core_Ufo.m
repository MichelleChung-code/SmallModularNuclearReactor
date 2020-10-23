function [outputArg1] = reactor_core_Ufo(fuel_elements_total_vol)
% Calculate the Overall Heat Transfer Coefficient

% CONSTANTS
Po = 250E6; % assume all the power is derived from heat
T_in = 250; 
T_out = 700;
D_fuel_elements = 6E-2; % assume transferred through surface area of spherical fuel elements (in m)
N_fuel_elements = 420000; %number of fuel elements in reactor core @ equilibrium state
% porosity = 0.39

%Ufo = Q/(A*deltaT)
A_particles = N_fuel_elements*(4 * pi * (D_fuel_elements/2)^2)
% SA = (1-porosity)*A_particles/fuel_elements_total_vol

Ufo = Po/(A_particles*(T_out-T_in));
outputArg1 = Ufo;
end

