function [outputArg1] = reactor_core_Ufo()
% Calculate the Overall Heat Transfer Coefficient

% CONSTANTS
Po = 250E6; % assume all the power is derived from heat
T_in = 250; 
T_out = 700;
D = 3; % assume transferred through cross sectional area

%Ufo = Q/(A*deltaT)
A = pi * (D/2)^2; % assume transferred through cross sectional area

Ufo = Po/(A*(T_out-T_in));
outputArg1 = Ufo;
end

