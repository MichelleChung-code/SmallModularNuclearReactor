clear, clear global

global DIST_PAR

DIST_PAR(1) = 1.5; % relative volatility
DIST_PAR(2) = 41; % total number of stages
DIST_PAR(3) =21; % feed stage 
DIST_PAR(4) = 1;  % feed flowrate 
DIST_PAR(5) = 0.5; % feed composition, light comp
DIST_PAR(6) = 1; % feed quality (1 = sat'd liqd,0 = sat'd vapor) 
DIST_PAR(7) = 2.706; % reflux
DIST_PAR(8) = 3.206; % reboiler vapor flowrate

ns = DIST_PAR(2);
disp(' ')

x0 = 0.99:-(0.98/(ns-1)):0.01; x0 = x0';
options = optimset('fsolve');
options.LargeScale = 'off';
[x] = fsolve('dist_ss',x0,options);
% x0 = x; save x0.mat x0;

DIST_PAR(7) = 2.706*1.01;
[xup] = fsolve('dist_ss',x0,options);

DIST_PAR(7) = 2.706*0.99;
[xdown] = fsolve('dist_ss',x0,options);

figure(3), clf, stages = 1:ns;
plot(stages,x,stages,xup,'r', stages,xdown,'g')
axis([1 ns 0 1])
legend('Base Case','+1% Reflux','-1% Reflux')
title('Steady-State Composition Profile'),
ylabel('Mole Fraction LK'), xlabel('Stage Number')

disp(' ')
save steady


  