% Random nice to have analytics
% Standalone script that can be run by itself
clc, close all, clear

fraction_of_rod_inserted(1).f1 = 0:.01:1;
relative_reactivities = cell2mat(arrayfun(@(x) x.f1 - (((2*pi)^-1)*sin(2*pi*x.f1)), fraction_of_rod_inserted,'UniformOutput',false));
fraction_of_rod_inserted = fraction_of_rod_inserted.f1;


figure(1), plot(fraction_of_rod_inserted, relative_reactivities), grid on 
title('Relative Reactivity Worth of Control Rods')
ylabel('\rho(x) / \rho(H)'), xlabel('x / H')
