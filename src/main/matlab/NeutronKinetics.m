classdef NeutronKinetics
    properties
        beta_ls_delayed_groups {mustBeNumeric}
        lambda_ls_delayed_groups {mustBeNumeric}
        coupling_coeffs_matrix {mustBeNumeric}
        beta {mustBeNumeric}
        lambda {mustBeNumeric}
        N {mustBeNumeric}
        alpha_fuel {mustBeNumeric}
        alpha_moderator {mustBeNumeric}
        alpha_reflector {mustBeNumeric}
    end
    
    methods
       function obj = NeutronKinetics(coupling_coeffs_matrix, N)
           % Constant class inputs that will not change 
           obj.beta_ls_delayed_groups = 10^-2*[.0256 .14 .13 .27 .086 .017]';
           obj.lambda = 7.66E-4;
           obj.beta = .67;
           obj.lambda_ls_delayed_groups = [.0256 .14 .13 .27 .086 .017]';
           
           % Dyanmic inputs that can change
           obj.coupling_coeffs_matrix = coupling_coeffs_matrix;
           obj.N = N; % number of nodes
           
           obj.alpha_fuel = alpha_fuel; 
           obj.alpha_moderator = alpha_moderator;
           obj.alpha_reflector = alpha_reflector;
       end
       
       function rho = reactivity(obj, Tc, Tr)
           % DOCUMENT WHAT EACH VALUE REPRESENT
           % Tc = Temperature of the fuel element
           % Tr = Temperature of the reflector
           % Tc0, tr0 = initial temperature of fuel elements and reflector 
           % rho_control_rods = reactivity introduced by the control rods
           
           rho_control_rods = 1.8e-3; % Need to find a valid number
           Tc0 = 1000; %Change number
           Tr0 = 1000; %Change number
           rho = rho_control_rods + (obj.alpha_fuel + obj.alpha_moderator)*(Tc - Tc0) + obj.alpha_reflector*(Tr -Tr0);
       
       end
       
       function dxdt = relative_neutron_flux_middle_node(obj, t, x, node)
           % DOCUMENT WHAT EACH VALUE REPRESENTS
           % x(1) = ni neutron flux of the ith node
           % x(2:7) = Concentrations of delayed groups for ith node
           
           n0 = 1 % ni-1 neutron flux of the ith-1 node
           n2 = 2 % ni+1 neutron flux of the ith+1 node
           rho_i = 1 % Reactivity of ith node 
           
           % Relative neutron flux for ith node
           dnidt_term1 = (rho_i - obj.beta - obj.coupling_coeffs_matrix(node,node))*(1/obj.lambda)*x(1);
           dnidt_term2 = (1/obj.lambda)*(obj.coupling_coeffs_matrix(node, node-1)*n0 + obj.coupling_coeff_matrix(node, node+1)*n2);
           dnidt_term3 = obj.sum_beta_concentration_over_lambda(x(2:7), obj.lambda);
           dxdt(1) = dnidt_term1 + dnidt_term2 + dnidt_term3;
           
           % Delayed Neutron Group Concentrations for Node i
           % TODO - see if there is a nicer way to do this - I hate how
           % this is done right now but solving the ODE45 inside the ODE45
           % was causing issues
           dxdt(2) = obj.lambda_ls_delayed_groups(1)*(x(1)-x(2)); 
           dxdt(3) = obj.lambda_ls_delayed_groups(2)*(x(1)-x(3));
           dxdt(4) = obj.lambda_ls_delayed_groups(3)*(x(1)-x(4));
           dxdt(5) = obj.lambda_ls_delayed_groups(4)*(x(1)-x(5));
           dxdt(6) = obj.lambda_ls_delayed_groups(5)*(x(1)-x(6));
           dxdt(7) = obj.lambda_ls_delayed_groups(6)*(x(1)-x(7));    
           
       end
       
       function dxdt = relative_neutron_flux_first_node(obj, t, x) 
           % DOCUMENT WHAT EACH VALUE REPRESENTS
           % x(1) = n1 neutron flux of the first node
           % x(2:7) = Concentrations of delayed groups for node 1
           % x(8) = nN neutron flux of the last node
           % x(9:14) = Concentrations of delayed groups for node N
           
           n2 = 1; % needs to be replaced with an equation 
           rho_1 = 1; % needs to be replaced with an equation
           rho_N = 1; % needs to be replaced with an equation
           n_N_1 = 1; % needs to be replaced with an equation
           
           % Relative neutron flux for node 1
           dn1dt_term1 = (rho_1 - obj.beta - obj.lambda_ls_delayed_groups(1,1))/obj.lambda*x(1);
           dn1dt_term2 = (1/obj.lambda) * obj.coupling_coeffs_matrix(1,2) * n2;
           dn1dt_term3 = obj.sum_beta_concentration_over_lambda(x(2:7), obj.lambda);
           dxdt(1) =  dn1dt_term1 + dn1dt_term2 + dn1dt_term3;
           
           % Delayed Neutron Group Concentrations for Node 1
           % TODO - see if there is a nicer way to do this - I hate how
           % this is done right now but solving the ODE45 inside the ODE45
           % was causing issues
           dxdt(2) = obj.lambda_ls_delayed_groups(1)*(x(1)-x(2)); 
           dxdt(3) = obj.lambda_ls_delayed_groups(2)*(x(1)-x(3));
           dxdt(4) = obj.lambda_ls_delayed_groups(3)*(x(1)-x(4));
           dxdt(5) = obj.lambda_ls_delayed_groups(4)*(x(1)-x(5));
           dxdt(6) = obj.lambda_ls_delayed_groups(5)*(x(1)-x(6));
           dxdt(7) = obj.lambda_ls_delayed_groups(6)*(x(1)-x(7));
           
           % Relative neutron flux for the last node (N)
           dnNdt_term1 = (rho_N - obj.beta - obj.coupling_coeffs_matrix(obj.N,obj.N))/obj.lambda*x(8);
           dnNdt_term2 = (1/obj.lambda)*obj.coupling_coeffs_matrix(obj.N,obj.N-1)*n_N_1;
           dnNdt_term3 = obj.sum_beta_concentration_over_lambda(x(9:14), obj.lambda);
           dxdt(8) =  dnNdt_term1 + dnNdt_term2 + dnNdt_term3;
           
           % Delayed Neutron Group Concentrations for Node N
           dxdt(9) = obj.lambda_ls_delayed_groups(1)*(x(8)-x(9)); 
           dxdt(10) = obj.lambda_ls_delayed_groups(2)*(x(8)-x(10));
           dxdt(11) = obj.lambda_ls_delayed_groups(3)*(x(8)-x(11));
           dxdt(12) = obj.lambda_ls_delayed_groups(4)*(x(8)-x(12));
           dxdt(13) = obj.lambda_ls_delayed_groups(5)*(x(8)-x(13));
           dxdt(14) = obj.lambda_ls_delayed_groups(6)*(x(8)-x(14));
           
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
       
       function [tout x] = solve_neutron_kinetics(obj, tspan, x0)
            % This function is going to solve the ODE 
            [tout, x] = ode45(@obj.relative_neutron_flux_first_node, tspan, x0);
            
       end

    end
    
end