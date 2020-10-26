classdef NeutronKinetics
    properties
        beta_ls_delayed_groups {mustBeNumeric}
        lambda_ls_delayed_groups {mustBeNumeric}
        coupling_coeffs_matrix {mustBeNumeric}
        beta {mustBeNumeric}
        lambda {mustBeNumeric}
    end
    
    methods
       function obj = NeutronKinetics(coupling_coeffs_matrix)
           % Constant class inputs that will not change 
           obj.beta_ls_delayed_groups = 10^-2*[.0256 .14 .13 .27 .086 .017]';
           obj.lambda = 7.66E-4;
           obj.beta = .67;
           obj.lambda_ls_delayed_groups = [.0256 .14 .13 .27 .086 .017]';
           
           % Dyanmic inputs that can change
           obj.coupling_coeffs_matrix = coupling_coeffs_matrix;
       end
       
       function dxdt = relative_neutron_flux_first_node(obj, t, x) 
           % x(1) = n1 neutron flux of the first node
           % x(2:7) = Concentrations of delayed groups for node 1
           % x(8) = nN neutron flux of the last node
           n2 = 1; % needs to be replaced with an equation 
           rho = 1; % needs to be replaced with an equation
           
           % Relative neutron flux for node 1
           dn1dt_term1 = (rho - obj.beta - obj.lambda_ls_delayed_groups(1,1))/obj.lambda_ls_delayed_groups(1)*x(1);
           dn1dt_term2 = (1/obj.lambda_ls_delayed_groups(1)) * obj.coupling_coeffs_matrix(1,2) * n2;
           dn1dt_term3 = obj.sum_beta_concentration_over_lambda(x(2:7), obj.lambda_ls_delayed_groups(1));
           dxdt(1) =  dn1dt_term1 + dn1dt_term2 + dn1dt_term3;
           
           % Delayed Neutron Group Concentrations for Node 1
           % TODO - see if there is a nicer way to do this - I hate how
           % this is done right now - MC
           dxdt(2) = obj.lambda_ls_delayed_groups(1)*(x(1)-x(2)); 
           dxdt(3) = obj.lambda_ls_delayed_groups(2)*(x(1)-x(3));
           dxdt(4) = obj.lambda_ls_delayed_groups(3)*(x(1)-x(4));
           dxdt(5) = obj.lambda_ls_delayed_groups(4)*(x(1)-x(5));
           dxdt(6) = obj.lambda_ls_delayed_groups(5)*(x(1)-x(6));
           dxdt(7) = obj.lambda_ls_delayed_groups(6)*(x(1)-x(7));
           
           % Relative neutron flux for the last node (N)
           
           
           dxdt = dxdt';
       end
       
       function sum_term3 = sum_beta_concentration_over_lambda(obj, C, lambda_node)
           % C is a vector containing 6 elements, the concentration of each
           % delayed neutron group, for a given node
           ls_elements_to_sum = zeros(1, length(C));
           for i=1:length(C)
               ls_elements_to_sum(i) = obj.beta_ls_delayed_groups(i) / lambda_node * C(i);
           end
       
           sum_term3 = sum(ls_elements_to_sum,'all');
  
       end
       
       function dCdt = relative_concentrations(obj, C, n)
           % concentration for one node
           % TODO
           % See if can figure out how to use this to minimize copy/pasting
           % code... issues in ODE45 within ODE45 function unfortunately -
           % MC
           dCdt(1) = obj.lambda_ls_delayed_groups(1)*(n-C(1));
           dCdt(2) = obj.lambda_ls_delayed_groups(2)*(n-C(2));
           dCdt(3) = obj.lambda_ls_delayed_groups(1)*(n-C(3));
           dCdt(4) = obj.lambda_ls_delayed_groups(2)*(n-C(4));
           dCdt(5) = obj.lambda_ls_delayed_groups(1)*(n-C(5));
           dCdt(6) = obj.lambda_ls_delayed_groups(2)*(n-C(6));
           
           dCdt = dCdt';
           
       end
       
       function [tout x] = solve_neutron_kinetics(obj, tspan, x0)
            % This function is going to solve the ODE 
            [tout, x] = ode45(@obj.relative_neutron_flux_first_node, tspan, x0);
            
       end

    end
    
end