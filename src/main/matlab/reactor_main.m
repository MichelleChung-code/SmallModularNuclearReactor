import pkg.*
% parameters of the delayed neutron group
half_life = 5.27;
beta = 0.006502; % effective delayed neutron fraction
lambda = 0.07669; % decay constant

reactor = SampleNuclearReactor(half_life, beta, lambda);

% for radioactive decay
[decay_time, decay_activity] = reactor.RadioActiveDecay();

% for point kinetic equation
[conc_time, concentration] = reactor.PointKineticEquation();






