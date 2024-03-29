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
        mass_flow_rate_helium {mustBeNumeric}
        Tout {mustBeNumeric}
        
        reactivity_step_size {mustBeNumeric} 
        reactivity_step_time {mustBeNumeric} 
        T_outlet_header_sp_step {mustBeNumeric} 
        T_outlet_header_sp_step_time {mustBeNumeric} 
        natural_reactivity {mustBeNumeric} 
        reactivity_profile {mustBeNumeric}

        T_outlet_header_ss {mustBeNumeric} 
        T_outlet_header_set_point {mustBeNumeric} 
        KC {mustBeNumeric} 
        TI {mustBeNumeric}
        tau {mustBeNumeric} 
        Kk {mustBeNumeric}
        controcost {mustBeNumeric}
        deadtime {mustBeNumeric}
        control_rod_insertion_ss {mustBeNumeric} 
        control_rod_min_insertion {mustBeNumeric} 
        control_rod_max_insertion {mustBeNumeric} 
        control_rod_reactivity {mustBeNumeric}
        
    end
    
    methods(Static)
        function vol = calc_volume(type, D, H)
           if nargin < 3
               H = 'Not Used'; % MATLAB does not support default arguments
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
       function obj = NeutronKinetics(coupling_coeffs_matrix, N, reactivity_step_size, reactivity_step_time, natural_reactivity,reactivity_profile,control_rod_reactivity, x0, control_rod_sp_step_time)
           % Dynamic inputs that can change
           obj.coupling_coeffs_matrix = coupling_coeffs_matrix;
           obj.N = N; % number of nodes
           obj.reactivity_step_size = reactivity_step_size; 
           obj.reactivity_step_time = reactivity_step_time; 
           obj.natural_reactivity = natural_reactivity;
           obj.reactivity_profile = reactivity_profile;
           obj.control_rod_reactivity = control_rod_reactivity;
           
           if obj.reactivity_step_size == 0
               disp("No step response inputted, running in Steady State Operation")
           end
           
           
           % VARIABLES DEFINED BELOW ARE CONSTANTS FOR THE
           % ENTIRE SYSTEM.
           % Dynamic inputs i.e. inputs that can change run to run should
           % be fed into the class object
           % Constant class inputs that will not change 
           obj.beta_ls_delayed_groups = 10^-2*[.0256 .14 .13 .27 .086 .017]';
           obj.lambda = 7.669E-2; 
           obj.beta = 10^-2*.67;
           obj.lambda_ls_delayed_groups = [.0124 .0305 .111 .301 1.14 3.01]'; 
           
           % Variable for reactivity
           obj.alpha_fuel = -4.36e-5; 
           obj.alpha_moderator = -0.94e-5;
           obj.alpha_reflector = 1.49e-5;
           obj.Tc0  = 0; % celsius - reference temperature  
           obj.Tr0 = 0; % celsius - reference temperature
           
           % Reactor Core Properties
           D_fuel_element = 6E-2; % m
           N_fuel_elements_in_core = 420000; % number of fuel elements
           D_reactor_core = 3; % m reactor core diameter
           H_reactor_core = 11; % m reactor core height
           reflector_thickness = 0.5; % m 
           
           %Constant for Thermal Hydraulics
           obj.volume = obj.calc_volume('cylinder', D_reactor_core, H_reactor_core); %m^3

           obj.porosity = 0.39;
           obj.density_fuel = 1797.169975; % kg/m^3 
           obj.volume_i = obj.volume/obj.N;
           
           obj.mass_flow_rate_helium = 145; % kg/s
           obj.Tin = 252.9; % Helium input temperature in Celsius 
           obj.Tout = 750; % Helium 
           Tavg_celsius = (obj.Tin+obj.Tout)/2;
           Tavg_Kelvin = Tavg_celsius + 273.15;
           Pin = 7000 / 1000; % MPa
           Pout = 6961.37 / 1000; % MPa
           Pavg = (Pin + Pout) / 2; % MPa
           
           % graphite thermal specific heat capacity using correlation from
           % http://aries.ucsd.edu/LIB/PROPS/PANOS/c.html
           % for 300 K <= T <= 1500 K
           
           obj.C_fuel = -474.0 + 4.9532*Tavg_Kelvin - 3.6093E-3*(Tavg_Kelvin^2) + 9.3068E-7*(Tavg_Kelvin^3); % j/kgK
           
%            Cp_helium from https://www-sciencedirect-com.ezproxy.lib.ucalgary.ca/science/article/pii/S0029549312000921
           obj.Cp_helium = 5195; % or 5246.5 J/kgK from Symmetry at 500 degrees Celsius 
           
           obj.P0 = obj.mass_flow_rate_helium * obj.Cp_helium * (obj.Tout - obj.Tin); %W
           obj.P0 = 7.602E8;%380371250;
           obj.P0_node = obj.P0/obj.N;  % Watts
           
           obj.density_reflector = 2230; % kg/m^3 (Density of graphite)file:Capstone_Group25_CHEMENGG/Reactor_Modelling/Reactor_Core_Modelling/Sources_Used/brochure-properties-and-characteristics-of-graphite-7329.pdf
           obj.volume_reflector = obj.calc_volume('cylinder', D_reactor_core + 2*(reflector_thickness), H_reactor_core); % m^3 
           
%            https://www-sciencedirect-com.ezproxy.lib.ucalgary.ca/science/article/pii/S0029549312000921
%            volume specific heat capacity for reflector graphite
           obj.C_reflector = 1.75*(0.645+3.14e-3*Tavg_celsius - 2.809e-6*(Tavg_celsius^2) +0.959e-9*(Tavg_celsius^3))*1000 %j/kgK
           obj.Ad = (1- obj.porosity)*(obj.calc_surface_area('sphere', D_fuel_element))*N_fuel_elements_in_core/obj.N; % m^2 SA of one sphere * number of spheres/10 sections
           obj.Kd = obj.P0/(obj.Ad*(obj.Tout-obj.Tin)); %116.9569; %W/(m^2*K) 
           obj.Kr = 70; % W/(m^2K) 
           obj.Ar = obj.calc_surface_area('cylinder_no_top', D_reactor_core, H_reactor_core)/obj.N *(1 - obj.porosity); % m^2 heat transfer area between fuel pile and reflector per node
           
           
           % Helium Properties 
           rho = 48.14*(Pavg/Tavg_Kelvin)*(1 + 0.4446*(Pavg/(Tavg_Kelvin^1.2)))^-1; % kg/m3
           mu = 3.674e-7*(Tavg_Kelvin^0.7);
           kinetic_visc = mu/rho; 
           
           SA_fuel_element = 4*pi*(D_fuel_element/2)^2;
           Re = ((obj.mass_flow_rate_helium * SA_fuel_element) * D_fuel_element)/ kinetic_visc;
           helium_heat_conductivity = 2.682e-3*(1+1.123e-3)*(Tavg_Kelvin^(0.71*(1-2e-4*Pavg)));% W/mK
           Pr = (mu * obj.Cp_helium) / helium_heat_conductivity;
           
           % heat transfer coefficient of the surface of spherical fuel elements
           Nu = (1.27*(Pr^(1/3)/obj.porosity^(1.18))*Re^(0.36)) + (0.033*(Pr^(1/2)/ obj.porosity^(1.07))*Re^(0.86));
           obj.K = helium_heat_conductivity * Nu / D_reactor_core; % W/(m^2K) heat transfer coefficient between fuel elements and helium
           obj.A = pi*(D_reactor_core/2)^2*(1- obj.porosity); % m^2 area of circle, cross-sectional area of the reactor core
           obj.Ku = 26494.5; % W/(m^2K) This number can be figure out by looking at the overall heat transfer between a flate plate and a flowing fluid
           obj.Au = obj.calc_surface_area('cylinder_no_top', .2,H_reactor_core)*30;%obj.calc_surface_area('cylinder_no_top', D_reactor_core + 2*(reflector_thickness), H_reactor_core); %m^2 heat transfer area between coolant in reflector and riser, SA of reflector using outer diameter
           
           obj.k = 0; %leakage ratio 

           % For the control rod control system
           obj.T_outlet_header_ss = obj.Tout; %Toh 
           obj.T_outlet_header_set_point = obj.Tout; 
           obj.T_outlet_header_sp_step = 25; % step the set point by 25 degrees
           
           % for no step change in set point, comment out the below line
           obj.T_outlet_header_sp_step = 0;
           
           obj.T_outlet_header_sp_step_time = control_rod_sp_step_time;
           obj.control_rod_insertion_ss = x0( obj.N*7+3*obj.N+6 + obj.N);  
           % Tuning parameters for PI control 
           obj.tau = 12; %time for Toh to reach .653 of the final output from a step change in reactivty
           obj.Kk = 650; %Change in Toh after a step change in reactivity
           obj.controcost = 11; % This is a guess
           obj.deadtime = 1;
           obj.KC = obj.tau/(obj.Kk*(obj.controcost +obj.deadtime));
           obj.TI = obj.tau;
           
           % control rod insertion must be between 0 and 11m 
           obj.control_rod_min_insertion = 0;  
           obj.control_rod_max_insertion = 11;
       end
       
       function rho = reactivity(obj,reactivity_profile, control_rod_x,Tc, Tr, t)
           % Args:
           % x = current control rod insertion position
           % Tc = Temperature of the fuel element
           % Tr = Temperature of the reflector
           % Tc0, tr0 = initial temperature of fuel elements and reflector 
           % rho_control_rods = reactivity introduced by the control rods
           
           % insertion position reference for future Ohki2014_Chapter_FuelBurnupAndReactivityControl
     
           rho_natural = obj.natural_reactivity*reactivity_profile;
           
           % STEP change in control rod natural reactivity
           if t >= obj.reactivity_step_time
               rho_natural = rho_natural + obj.reactivity_step_size;
           end
           
           rho_control_rods = -control_rod_x*obj.control_rod_reactivity; 
           rho = rho_natural +rho_control_rods + (obj.alpha_fuel + obj.alpha_moderator)*(Tc - obj.Tc0) + obj.alpha_reflector*(Tr - obj.Tr0);

       end
       function dxdt = relative_neutron_flux(obj, t, x) 
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
           integ = x(length(x)); % for PI controller

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control_rod_insertion is the Manipulated variable
% Toh is the controlled variable
% PI control
           
           global reactivity_control_rod_results
           global time_plot
            
           % step change in the set point 
           if t >= obj.T_outlet_header_sp_step_time
               obj.T_outlet_header_set_point = obj.T_outlet_header_ss + obj.T_outlet_header_sp_step;
           end

           error = x(Toh) - obj.T_outlet_header_set_point; 
           control_rod_insertion = obj.control_rod_insertion_ss + obj.KC*(error+1/obj.TI*integ);

           % Clamp the manipulated variable
           control_rod_insertion = min(control_rod_insertion, obj.control_rod_max_insertion);
           control_rod_insertion = max(control_rod_insertion, obj.control_rod_min_insertion);
           
           dxdt(length(x)) = error; % integ
           dxdt(control_rod_position) = control_rod_insertion;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           H_reactor_node = 11/obj.N;
           controller_node_number = control_rod_insertion/H_reactor_node;
           
           control_rod_array = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
           
           %Calculation for percent control rod insertion per node%
           for i = 1: obj.N
               if controller_node_number >= 1
                   control_rod_array(i) = 1;
                   controller_node_number = controller_node_number-1;
               elseif controller_node_number > 0
                   control_rod_array(i) = controller_node_number;
                   controller_node_number = controller_node_number-1;
               end
           end
           control_rod_array = control_rod_array.*obj.reactivity_profile;

           
           % Relative neutron flux for node 1
           rho_1 = obj.reactivity(obj.reactivity_profile(1),control_rod_array(1),x(cores),x(Tr), t);
           dxdt(rho_index) = rho_1;
           dn1dt_term1 = (rho_1 - obj.beta - obj.coupling_coeffs_matrix(1,1))/obj.lambda*x(1);
           dn1dt_term2 = (1/obj.lambda) * obj.coupling_coeffs_matrix(1,2) * x(2);
           dn1dt_term3 = obj.sum_beta_concentration_over_lambda(x(obj.N+1:obj.N+6), obj.lambda);
           dxdt(1) =  dn1dt_term1 + dn1dt_term2 + dn1dt_term3;
           
           intermediate_rho_ls = [];
           % Relative neutron flux for ith to N-1 nodes
           var = obj.N+7;
           for i = 2:(obj.N-1)
               rho_i = obj.reactivity(obj.reactivity_profile(i),control_rod_array(i), x(cores+i-1),x(Tr), t); 
               dxdt(rho_index + i -1) = rho_i;
               dnidt_term1 = (rho_i - obj.beta - obj.coupling_coeffs_matrix(i,i))*(1/obj.lambda)*x(i);
               dnidt_term2 = (1/obj.lambda)*(obj.coupling_coeffs_matrix(i, i-1)*x(i-1) + obj.coupling_coeffs_matrix(i, i+1)*x(i+1));
               dnidt_term3 = obj.sum_beta_concentration_over_lambda(x(var:var+5), obj.lambda);
               dxdt(i) = dnidt_term1 + dnidt_term2 + dnidt_term3;
               var = var+6; 
               intermediate_rho_ls = [intermediate_rho_ls, rho_i];
           end
           
           % Relative neutron flux for node N
           rho_N = obj.reactivity(obj.reactivity_profile(obj.N),control_rod_array(obj.N), x(cores+obj.N-1),x(Tr), t); 
           dxdt(rho_index + obj.N -1) = rho_N; 
           dnNdt_term1 = (rho_N - obj.beta - obj.coupling_coeffs_matrix(obj.N,obj.N))/obj.lambda*x(obj.N);
           dnNdt_term2 = (1/obj.lambda)*obj.coupling_coeffs_matrix(obj.N,obj.N-1)*x(obj.N-1);
           dnNdt_term3 = obj.sum_beta_concentration_over_lambda(x(obj.N*7-5:obj.N*7), obj.lambda);
           dxdt(obj.N) =  dnNdt_term1 + dnNdt_term2 + dnNdt_term3;
          
           
           % Output reactivities and control rod insertion into a global
           % variable
           reactivity_control_rod_results = [reactivity_control_rod_results; [rho_1 intermediate_rho_ls rho_N control_rod_insertion]];
           time_plot = [time_plot t];
           
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toh is the controlled variable 743.8 = set point 
           dxdt(Toh) = (dTohdt_term1 + dTohdt_term2)/Woh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
           
%            Initiate global variable for plotting reactivity and control
%            rod insertion 

            global reactivity_control_rod_results
            global time_plot

           
            % Function called to simultaneously solve system of
            % differential equations
            
            tic
            tstep = .1;
            tspan_fix = tspan(1):tstep:tspan(2);
            reactivity_control_rod_results = [];
            time_plot = [];
            [tout, x] = ode23s(@obj.relative_neutron_flux, tspan_fix, x0);
            
            toc
            
            % POST-PROCESSING
            % Divide nodal reactivity by tout, as it was included in the
            % ODE solving to extract reactivity values and integrated over
            % time unneccessarily when wanting to extract for plotting
            % purposes
            
            [num_row_x, num_col_x] = size(x);
            global num_col_x_without_PI
            num_col_x_without_PI = num_col_x -1;
            reactivity_final_col_num = num_col_x_without_PI - 1;
            x(:, reactivity_final_col_num-(obj.N-1):reactivity_final_col_num) = rdivide(x(:, reactivity_final_col_num-(obj.N-1):reactivity_final_col_num), tout);
            
            x(:, reactivity_final_col_num + 1) = rdivide(x(:, reactivity_final_col_num + 1), tout);
            % calculate power output per node
            x(:, num_col_x_without_PI + 1:num_col_x_without_PI + obj.N) = x(:, 1:obj.N) * (1/obj.N) * (obj.P0*10^-6);
            
            % normalize the concentrations with the SS of the first node
%             [g1_fact,g2_fact,g3_fact,g4_fact,g5_fact,g6_fact] = subsref(num2cell(x(length(tout), obj.N + 1:obj.N + 6)),struct('type',{'{}'},'subs',{{1:6}}));
%             LS_normalize_facts = [g1_fact,g2_fact,g3_fact,g4_fact,g5_fact,g6_fact];
%             counter = obj.N;
%             for j = 1:obj.N
%               for i = 1:6
%                   x(:, i+counter) = x(:, i+counter)/LS_normalize_facts(i);
%               end
%               counter = counter + 6;
%             end
       end

    end
    
end
