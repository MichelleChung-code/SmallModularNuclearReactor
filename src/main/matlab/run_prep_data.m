close all, clear
clc
% Just to call functions

[fuelMass, fuelVolume] = fuel_mass_vol() % in kg, in m^3
Ufo = reactor_core_Ufo() % in W/(m^2*K)

% Returns: 
% Fuel Mass: 8.5367e+04 kg
% Fuel Volume: 47.5009 m^3
% Ufo: 7.8595e+04 W/(m^2*K)


