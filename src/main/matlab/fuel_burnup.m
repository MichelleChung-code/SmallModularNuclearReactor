function natural_reactivity = fuel_burnup(avg_fuel_burnup, num_nodes, natural_reactivity)
% All at equilibrum rate
% Discharge burnup val = 80000 MWd/tU

% https://www-sciencedirect-com.ezproxy.lib.ucalgary.ca/science/article/pii/S0029549321000157
% 80% of discharged fuel pebbles were put back into core and other 20%
% removed and replaced by fresh fuel element

    % per node fuel burnup change
    burn_up_change = avg_fuel_burnup/ 1000 / (num_nodes - 1);
    % maybe weight by loading amount instead
    
    

    % array of impact per node
    
    
    % initiate array of zeros 
    fuel_burnup_val = zeros(1, 10);
    fuel_burnup_res = zeros(1, 10);
    
    for i=1 : length(fuel_burnup_res)
        if i > 1
            fuel_burnup_val(i) = fuel_burnup_val(i-1) + burn_up_change;
        end
        x = fuel_burnup_val(i);
        fuel_burnup_res(i) = 3E-8*x^4 - 4E-6*x^3 + 0.0001*x^2 - 0.0045*x + 1.1148; % keff
        
        fuel_burnup_res(i) = (fuel_burnup_res(i)-1) / fuel_burnup_res(i); %rho
        
        % decrease due to enrichment differences between literature source
        % and SMR
        % burnup profile literature source:
        % https://www.sciencedirect.com/science/article/pii/S1738573319303043
        fuel_burnup_res(i) =  fuel_burnup_res(i) * (1 - (17-8.77)/100);
        
        % added impact
        natural_reactivity(i) =  fuel_burnup_res(i) + natural_reactivity(i);
        
    end
    
    
end