classdef SampleNuclearReactor
    properties
        half_life {mustBeNumeric} % half_life of the radioactive source
        base_run_path
        simulink_run_path
        beta {mustBeNumeric}
        lambda {mustBeNumeric}
    end
    
    methods
       function obj = SampleNuclearReactor(half_life, beta, lambda)
           obj.half_life = half_life;
           obj.base_run_path = fileparts(pwd); 
           obj.simulink_run_path = strcat(obj.base_run_path, '\simulink\');
           obj.beta = beta; % effective delayed neutron fraction
           obj.lambda = lambda; % decay constant of the delayed neutron group
       end
       
       function [decay_time, decay_activity] = RadioActiveDecay(obj)
           decay_run_path = strcat(obj.simulink_run_path, 'radioactive_decay.slx');
           sim(decay_run_path); 
           decay_time = tout; % decay_time in years
           decay_activity = activity; % Decay activity
       end
       
       function [conc_time, concentration] = PointKineticEquation(obj)
           % point kinetic equation based on treatise by Hetrick
           point_kin_run_path = strcat(obj.simulink_run_path, 'point_kinetic_equation.slx');
           sim(point_kin_run_path); 
           conc_time = tout;
           concentration = conc;
           
       end 
   end
       
end