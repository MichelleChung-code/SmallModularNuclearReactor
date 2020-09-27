classdef SampleNuclearReactor
    properties
        half_life {mustBeNumeric} % half_life of the radioactive source
        base_run_path
    end
    
    methods
       function obj = SampleNuclearReactor(half_life)
           obj.half_life = half_life;
           obj.base_run_path = fileparts(pwd); 
       end
       
       function RadioActiveDecay(obj)
           simulink_run_path = strcat(string(obj.base_run_path), '\simulink\radioactive_decay.slx');
           sim(simulink_run_path)
           
       end
   end
       
end