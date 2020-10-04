close all, clear
clc
% Just to call functions
d_fuel_kernel = 2*0.25e-3; % in m
ls_coating_layer_t = [0, 0.09e-3, 0.04e-3, 0.035e-3, 0.04e-3]; % in m
ls_coating_layer_rho = [10.4, 1.1, 1.9, 3.18, 1.9];
ls_coating_layer_rho = ls_coating_layer_rho*100^3; % in g/(m^3)


fuelMass = fuel_mass(d_fuel_kernel, ls_coating_layer_t, ls_coating_layer_t) % in kg

% Returns: 8.9864e+04

