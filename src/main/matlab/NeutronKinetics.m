classdef NeutronKinetics
    properties
        beta_ls_delayed_groups {mustBeVector} 
        lambda_ls_delayed_groups {mustBeVector}
        coupling_coeffs_matrix {mustBeVector}
        beta {mustBeNumeric}
        lambda {mustBeNumeric}
        tspan {mustBeVector}
        x0 {mustBeVector}
    end
    
    methods
       function obj = NeutronKinetics(coupling_coeffs_matrix, tspan, x0)
           % Constant class inputs that will not change 
           obj.beta_ls_delayed_groups = 10^-2*[.0256 .14 .13 .27 .086 .017]';
           obj.lambda = 7.66E-4;
           obj.beta = .67;
           obj.lambda_ls_delayed_groups = [.0256 .14 .13 .27 .086 .017]';
           obj.tspan = tspan;
           obj.x0 = x0 
           
           % Dyanmic inputs that can change
           obj.coupling_coeffs_matrix = coupling_coeffs_matrix;
       end
       
       function dxdt = relative_neutron_flux_first_node(obj) 
           % x(1) = n1
           % X(2) = nN (the last node) 
           % 
           dn1dt_term1 = (rho - obj.beta - obj.lambda_ls_delayed_groups(1,1))/lambda_ls_delayed_groups(1)*n1
           dn1dt_term2 = (1/obj.lambda_ls_delayed_groups(1) * obj.coupling_coeffs_matrix(1,2) * 
           
       end
       
       function [tout x] = solve_neutron_kinetics(obj)
            % This function is going to solve the ODE 
       end

    end
    
end