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
        volume {mustBeNumeric}
        P0 {mustBeNumeric}
        P0_node {mustBeNumeric}
        porosity {mustBeNumeric}
        density_fuel {mustBeNumeric}
        volume_i {mustBeNumeric}
        C_fuel {mustBeNumeric}
        Cp_helium {mustBeNumeric}
        density_reflector {mustBeNumeric}
        volume_reflector {mustBeNumeric}
        C_reflector {mustBeNumeric}
        Tc0 {mustBeNumeric}
        Tr0 {mustBeNumeric}
        
        Kd {mustBeNumeric}
        Ad {mustBeNumeric}
        Kr {mustBeNumeric}
        Ar {mustBeNumeric}
        K {mustBeNumeric}
        A {mustBeNumeric}
        Ku {mustBeNumeric}
        Au {mustBeNumeric}
        k {mustBeNumeric}
        Tin {mustBeNumeric}
        control_rod_length {mustBeNumeric}
        mass_flow_rate_helium {mustBeNumeric}
        Tout {mustBeNumeric}
        
    end
    
    methods(Static)
        function vol = calc_volume(type, D, H)
           if nargin < 3
               H = 'Not Used'; % apparently MATLAB doesn't support default arguments XD
           end
           if strcmp('sphere', type)
              vol =  (4/3)*pi*(D/2)^2;
           elseif strcmp('cylinder', type)
              vol = pi*(D/2)^2*H;
           else
              error("Input type for volume calculation is not implemented")
           end
        end
        
        function SA = calc_surface_area(type, D, H)
           if nargin < 3
              H = 'Not Used'; 
           end
           if strcmp('sphere', type)
              SA =  4*pi*(D/2)^2;
           elseif strcmp('cylinder', type)
              SA = (2*pi*(D/2)*H) + 2*pi*(D/2)^2;
           elseif strcmp('cylinder_no_top', type)
              SA = (2*pi*(D/2)*H);
           else
              error("Input type for SA calculation is not implemented")
           end
        end
        
    end
    
    methods
       function obj = NeutronKinetics(coupling_coeffs_matrix, N)
           % Dyanmic inputs that can change
           obj.coupling_coeffs_matrix = coupling_coeffs_matrix;
           obj.N = N; % number of nodes
           
           % ALL VARIABLES THAT ARE DEFINED BELOW ARE CONSTANTS FOR THE
           % ENTIRE SYSTEM.
           % Dynamic inputs i.e. inputs that can change run to run should
           % be fed into the class object
           % Constant class inputs that will not change 
           obj.beta_ls_delayed_groups = 10^-2*[.0256 .14 .13 .27 .086 .017]';
           obj.lambda = 7.669E-2; % uses initial paper's values
           obj.beta = 10^-2*.67;
           obj.lambda_ls_delayed_groups = [.0124 .0305 .111 .301 1.14 3.01]'; % uses initial paper's values
           
           %Variable for reactivity
           obj.alpha_fuel = -4.36e-5; 
           obj.alpha_moderator = -0.94e-5;
           obj.alpha_reflector = 1.49e-5;
           obj.Tc0  = 0; %celsius CAN MODIFY
           obj.Tr0 = 0; %celsius CAN MODIFY
           
           % Reactor Core Properties
           D_fuel_element = 6E-2; % m
           N_fuel_elements_in_core = 420000; % number of fuel elements
           D_reactor_core = 3; % m reactor core diameter
           H_reactor_core = 11; % m reactor core height
           reflector_thickness = 0.5; % m we need to actually find this value CAN MODIFY
           
           %Constant for Thermal Hydraulics
           obj.volume = obj.calc_volume('cylinder', D_reactor_core, H_reactor_core); %m^3

           obj.porosity = 0.39;
           obj.density_fuel = 1797.169975; % kg/m^3 From stream table of our PFD
           obj.volume_i = obj.volume/obj.N;
           obj.C_fuel = 1621.45; % j/kgK
           obj.Cp_helium = 5.19E3; % j/kgK
           obj.mass_flow_rate_helium = 145; % kg/s
           obj.Tin = 250; % Helium input temperature in Celsius 
           obj.Tout = 700; % Helium 
           obj.P0 = obj.mass_flow_rate_helium * obj.Cp_helium * (obj.Tout - obj.Tin); %W
           obj.P0_node = obj.P0/obj.N;  % Watts
           
           
           %Values that need to be defined for now these are all random
           %number I created
           obj.density_reflector = 1760; % kg/m^3 CAN MODIFY (Density of graphite)file:///C:/Users/marce/University%20of%20Calgary/Michelle%20Chung%20-%20Capstone_Group25_CHEMENGG/Reactor_Modelling/Reactor_Core_Modelling/Sources_Used/brochure-properties-and-characteristics-of-graphite-7329.pdf
           obj.volume_reflector = obj.calc_volume('cylinder', D_reactor_core + 2*(reflector_thickness), H_reactor_core); % m^3 
           obj.C_reflector = 1735.5; % j/kgK CAN MODIFY (value of specific heat capacity of graphite at 750C)    file:///C:/Users/marce/University%20of%20Calgary/Michelle%20Chung%20-%20Capstone_Group25_CHEMENGG/Reactor_Modelling/Reactor_Core_Modelling/Sources_Used/brochure-properties-and-characteristics-of-graphite-7329.pdf
           obj.Ad = (1- obj.porosity)*(obj.calc_surface_area('sphere', D_fuel_element))*N_fuel_elements_in_core/obj.N; % m^2 SA of one sphere * number of spheres/10 sections
           obj.Kd = obj.P0/(obj.Ad*(obj.Tout-obj.Tin)); %116.9569; %W/(m^2*K) include porosity
           obj.Kr = 60; % W/(m^2K) CAN MODIFY 
           obj.Ar = obj.calc_surface_area('cylinder_no_top', D_reactor_core, H_reactor_core)/obj.N *(1 - obj.porosity); % m^2 heat transfer area between fuel pile and reflector per node
           obj.K = 60; % W/(m^2K) CAN MODIFY
           obj.A = pi*(D_reactor_core/2)^2*(1- obj.porosity); % m^2 area of circle, cross-sectional area of the reactor core
           obj.Ku = 112; % W/(m^2K) look into what the material is, fluid and wall CAN MODIFY This number can be figure out by looking at the overall heat transfer between a flate plate and a flowing fluid!!!!!! Should be ho 
           obj.Au = obj.calc_surface_area('cylinder_no_top', D_reactor_core + 2*(reflector_thickness), H_reactor_core); %m^2 heat transfer area between coolant in reflector and riser, SA of reflector using outer diameter
           
           obj.k = 0; %leakage ratio CAN MODIFY
           obj.control_rod_length = 4; % I made this up.  The HTR-10 value was 2.2m.  We need to find the HTR-PM value
       end
       
       function rho = reactivity(obj, control_rod_x,Tc, Tr)
           % Args:
           % x = current control rod insertion position
           % Tc = Temperature of the fuel element
           % Tr = Temperature of the reflector
           
           % Tc0, tr0 = initial temperature of fuel elements and reflector 
           % rho_control_rods = reactivity introduced by the control rods
           
           % control rod reactivity is defined as a sin wave based o1n the
           % insertion position Ohki2014_Chapter_FuelBurnupAndReactivityControl
           
           control_rod_fraction_inserted = control_rod_x / obj.control_rod_length; 
           rho_control_rods_H = 5.25E-2; % from "Capstone_Group25_CHEMENGG\Reactor_Modelling\2006-design-aspect-of-the-chinese-modular-high-temperature-gas-cooled-reactor-htr-pm_zhang.pdf" control rod worth
           rho_control_rods = rho_control_rods_H * (control_rod_fraction_inserted - (1/(2*pi))*sin(2*pi*control_rod_fraction_inserted)); 
           rho_control_rods = 0.028; % CAN CHANGE TO SEE IF BETTER GRAPH
           rho = rho_control_rods + (obj.alpha_fuel + obj.alpha_moderator)*(Tc - obj.Tc0) + obj.alpha_reflector*(Tr - obj.Tr0);

       end
       function dxdt = relative_neutron_flux(obj, ~, x) 
           % Function solving simulataneous differential equations
           % Args:
           % ~ = timespan to solve over
           % x = array containing intiial conditions for each differential
           % equation
           
           % DOCUMENT WHAT EACH X VALUE REPRESENTS
           % x(1:10) = ni neutron flux of the nodes
           % x(11:70) = Concentrations of delayed groups for nodes
           % x(71:80) = ith Temperature of Fuel element nodes Tci
           % x(81:90) = ith Temperature of Helium nodes Tdi
           Tr = obj.N*7+2*obj.N+1; % x(91) = Temperature of Reflector Tr
           Tu = obj.N*7+2*obj.N+2; % x(92) = Average Temperature in riser Tu
           Tlh = obj.N*7+2*obj.N+3; % x(93) = Temperature of lower helium header Tlh
           Toh = obj.N*7+2*obj.N+4; % x(94) = Temperature of outlet header Toh
           % x(95:104) = mass flowrate of helium at ith nodes h_mass
           Wu = obj.N*7+3*obj.N+5; % x(105) = mass flowrate of riser Wu
           % x(106:115) = reactivity values, these are not differential
           % equations and are just included for plotting purposes
           rho_index = Wu + 1; % x(106) starting index to go to the reactivity values
           control_rod_position = rho_index + obj.N; % x(116) starting index to go to the control rod insertion position value
           dxdt = zeros(1, control_rod_position); % predefine size of result array for performance
           Wlh = obj.mass_flow_rate_helium; % kg/s = mass flowrate for lower helium header Wlh
           cores = obj.N*7+1; % This is the index number for core (fuel element) temperature
           downs = obj.N*7+obj.N+1; % This is the index number for downcomer (helium) temperature
           hmass = obj.N*7+2*obj.N+5; % index number for masses

           % control rod position - intech-the_theoretical_simulation_of_a_model_by_simulink_for_surveying_the_work_and_dynamical_stability_of_nuclear_reactors_cores (1)
           control_rod_const_coeff = 0.1; % I made this up, this is the constant coefficient CAN MODIFY
           Ko = 0.5; % I made this up, this is the initial value of Keff CAN MODIFY
           Ksp = 0.5; % I made this up, this is supposed to be the secondary value of Keff in the recent position of control rod... Need to figure out how to access previous results in matlab ODE solver to get this
           velocity_control_rod = 10; % Units of mm/s I made this up too, this is the control rod velocity CAN MODIFY
           F = x(control_rod_position)*obj.control_rod_length*control_rod_const_coeff + Ko;
           dxdt(control_rod_position) = velocity_control_rod*sign(F - Ksp);
        
           
           % Relative neutron flux for node 1
           rho_1 = obj.reactivity(x(control_rod_position),x(cores),x(Tr));
           dxdt(rho_index) = rho_1;
           dn1dt_term1 = (rho_1 - obj.beta - obj.coupling_coeffs_matrix(1,1))/obj.lambda*x(1);
           dn1dt_term2 = (1/obj.lambda) * obj.coupling_coeffs_matrix(1,2) * x(2);
           dn1dt_term3 = obj.sum_beta_concentration_over_lambda(x(obj.N+1:obj.N+6), obj.lambda);
           dxdt(1) =  dn1dt_term1 + dn1dt_term2 + dn1dt_term3;
           
           % Relative neutron flux for ith to N-1 nodes
           var = obj.N+7;
           for i = 2:(obj.N-1)
               rho_i = obj.reactivity(x(control_rod_position), x(cores+i-1),x(Tr)); 
               dxdt(rho_index + i -1) = rho_i;
               dnidt_term1 = (rho_i - obj.beta - obj.coupling_coeffs_matrix(i,i))*(1/obj.lambda)*x(i);
               dnidt_term2 = (1/obj.lambda)*(obj.coupling_coeffs_matrix(i, i-1)*x(i-1) + obj.coupling_coeffs_matrix(i, i+1)*x(i+1));
               dnidt_term3 = obj.sum_beta_concentration_over_lambda(x(var:var+5), obj.lambda);
               dxdt(i) = dnidt_term1 + dnidt_term2 + dnidt_term3;
               var = var+6; 
           end
           
           % Relative neutron flux for node N
           rho_N = obj.reactivity(x(control_rod_position), x(cores+obj.N-1),x(Tr)); 
           dxdt(rho_index + obj.N -1) = rho_N; 
           dnNdt_term1 = (rho_N - obj.beta - obj.coupling_coeffs_matrix(obj.N,obj.N))/obj.lambda*x(obj.N);
           dnNdt_term2 = (1/obj.lambda)*obj.coupling_coeffs_matrix(obj.N,obj.N-1)*x(obj.N-1);
           dnNdt_term3 = obj.sum_beta_concentration_over_lambda(x(obj.N*7-5:obj.N*7), obj.lambda);
           dxdt(obj.N) =  dnNdt_term1 + dnNdt_term2 + dnNdt_term3;
          
           % Delayed Neutron Group Concentrations for Nodes
           var = obj.N;
           for i = 1:obj.N
              for j = 1:6
                  dxdt(var+j) = obj.lambda_ls_delayed_groups(j)*(x(i)-x(var+j));
              end
              var = var + 6;
           end
           

           % THERMAL HYDROLICS         
           % Fuel Pile Temperature for 1st node
           dTc1dt_term1 = obj.P0_node*x(1);
           dTc1dt_term2 = - obj.Kd*obj.Ad*(x(cores)-x(downs));
           dTc1dt_term3 = - obj.Kr*obj.Ar*(x(cores)-x(Tr));
           dTc1dt_term4 = -obj.K*obj.A*(x(cores)-x(cores+1));
           dxdt(cores) = (dTc1dt_term1 +dTc1dt_term2 +dTc1dt_term3+ dTc1dt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);
           
           % Fuel Pile Temperature for ith to N-1 node
           for i = 1:obj.N-2
               dTcidt_term1 = obj.P0_node*x(i+1);
               dTcidt_term2 = - obj.Kd*obj.Ad*(x(cores+i)-x(downs+i));
               dTcidt_term3 = - obj.Kr*obj.Ar*(x(cores+i)-x(Tr));
               dTcidt_term4 = -obj.K*obj.A*(x(cores+i)-x(cores+i+1)) + obj.K*obj.A*(x(cores+i-1)- x(cores+i));
               dxdt(cores+i) = (dTcidt_term1 +dTcidt_term2 +dTcidt_term3 + dTcidt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);               
           end
           
           % Fuel Pile Temperature for Nth node
           dTcNdt_term1 = obj.P0_node*x(obj.N);
           dTcNdt_term2 = - obj.Kd*obj.Ad*(x(cores + obj.N-1)-x(downs + obj.N-1));
           dTcNdt_term3 = - obj.Kr*obj.Ar*(x(cores + obj.N-1)-x(Tr));
           dTcNdt_term4 = obj.K*obj.A*(x(cores + obj.N-1 - 1)-x(cores+ + obj.N-1));
           dxdt(cores + obj.N-1) = (dTcNdt_term1 +dTcNdt_term2 +dTcNdt_term3+ dTcNdt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);           
           
           % Reflector Temperature           
           dTrdt_term1 = 0;
           for i = 0:obj.N -1
              dTrdt_term1 = dTrdt_term1 + obj.Kr*obj.Ar*(x(cores +i)-x(Tr)); 
           end
           dTrdt_term2 = -obj.Ku*obj.Au*(x(Tr)-x(Tu));
           dxdt(Tr) = (dTrdt_term1 + dTrdt_term2)/(obj.density_reflector*obj.volume_reflector*obj.C_reflector);
           
           % Downcomer Temperatrues
           %hmass node 1
           dxdt(hmass) = x(Wu) - x(hmass);
           %density node 2:N
           for i = 1:obj.N-1
               dxdt(hmass+i) = x(hmass+i-1) - x(hmass+i);
           end
           
           %Td node 1
           dTd1dt_term1 = obj.Kd*obj.Ad*(x(cores)- x(downs));
           dTd1dt_term2 = obj.Cp_helium*x(Wu)*(x(Tu)- x(downs));
           dxdt(downs) = (dTd1dt_term1 + dTd1dt_term2)/(obj.porosity*x(hmass)*obj.Cp_helium);
           
           %Td nodes 2:N
           for i = 1:obj.N-1
               dTdidt_term1 = obj.Kd*obj.Ad*(x(cores+i)-x(downs+i));
               dTdidt_term2 = obj.Cp_helium*x(hmass+i-1)*(x(downs+i-1)-x(downs+i));
               dxdt(downs+i)= (dTdidt_term1 + dTdidt_term2)/(obj.porosity*x(hmass+i)*obj.Cp_helium);
           end
           
           %Riser Temperature
           %hmass of riser Wu
           dxdt(Wu) = (1-obj.k)*Wlh - x(Wu);
           
           %Tu 
           dTudt_term1 = obj.Cp_helium*(1-obj.k)*Wlh*(x(Tlh) - x(Tu));
           dTudt_term2 = obj.Ku*obj.Au*(x(Tr) - x(Tu));
           dxdt(Tu) = (dTudt_term1 + dTudt_term2)/(x(Wu)*obj.Cp_helium);
           
           %Lower Header Tlh
           dxdt(Tlh) = (obj.Tin - x(Tlh));
           
           %Outer Header Toh
           Woh = obj.k*Wlh + x(hmass+obj.N-1);
           
           dTohdt_term1 = obj.k*Wlh*(x(Tlh)-x(Toh));
           dTohdt_term2 = x(hmass+obj.N-1)*(x(downs+obj.N-1)-x(Toh));
           dxdt(Toh) = (dTohdt_term1 + dTohdt_term2)/Woh;
          
           dxdt = dxdt';
       end
  
       
       function result = sum_beta_concentration_over_lambda(obj, C, lambda_node)
           % C is a vector containing 6 elements, the concentration of each
           % delayed neutron group, for a given node
           ls_elements_to_sum = zeros(1, length(C));
           for i=1:length(C)
               ls_elements_to_sum(i) = obj.beta_ls_delayed_groups(i) / lambda_node * C(i);
           end
           result = sum(ls_elements_to_sum,'all');
       end
       
       function [tout, x] = solve_neutron_kinetics(obj, tspan, x0)
            % Function called to simultaneously solve system of
            % differential equations
            
            tic
            tstep = .1;
            tspan_fix = tspan(1):tstep:tspan(2);
            [tout, x] = ode23s(@obj.relative_neutron_flux, tspan_fix, x0);
            
            toc
            
            % POST-PROCESSING
            % Don't hardcode this 
            % Divide nodal reactivity by tout, as it was included in the
            % ODE solving to extract reactivity values
            x(:, 106:115) = rdivide(x(:, 106:115), tout);
            
            % calculate power output
            % Per node
            x(:, 117:126) = x(:, 1:10) * (1/obj.N) * (obj.P0*10^-6);
            
            % normalize the concentrations with the SS of the first node
            [g1_fact,g2_fact,g3_fact,g4_fact,g5_fact,g6_fact] = subsref(num2cell(x(length(tout), 11:16)),struct('type',{'{}'},'subs',{{1:6}}));
            LS_normalize_facts = [g1_fact,g2_fact,g3_fact,g4_fact,g5_fact,g6_fact];
            counter = obj.N;
            for j = 1:obj.N
                for i = 1:6
                    x(:, i+counter) = x(:, i+counter)/LS_normalize_facts(i);
                end
                counter = counter + 6;
            end
       end

    end
    
end
