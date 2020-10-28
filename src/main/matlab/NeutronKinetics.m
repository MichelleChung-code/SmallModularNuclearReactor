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
           
           %Variable for reactivity
           obj.alpha_fuel = -4.36e-5; 
           obj.alpha_moderator = -0.94e-5;
           obj.alpha_reflector = 1.49e-5;
           obj.Tc0  = 1000; %celsius
           obj.Tr0 = 475; %celsius 
           
           %Constant for Thermal Hydraulics
           obj.volume = 77.75441818; %m^3
           obj.P0 = 250E6;  %Watts
           obj.porosity = 0.39;
           obj.density_fuel = 1797.169975; % kg/m^3 From strem table of our PFD
           obj.volume_i = obj.volume/obj.N;
           obj.C_fuel = 1621.45; %j/kgK
           obj.Cp_helium = 5.19E-3; % j/kgK
           %The following I just made up
           obj.density_reflector = 2000; % kg/m^3
           obj.volume_reflector = 112.31193; % m^3 assuming that reflector all around and thickness of .5m
           obj.C_reflector = 1621.45; % j/kgK            
       end
       
       function rho = reactivity(obj, Tc, Tr)
           % DOCUMENT WHAT EACH VALUE REPRESENT
           % Tc = Temperature of the fuel element
           % Tr = Temperature of the reflector
           % Tc0, tr0 = initial temperature of fuel elements and reflector 
           % rho_control_rods = reactivity introduced by the control rods
           
           rho_control_rods = 1.8e-3; % Need to find a valid number
           rho = rho_control_rods + (obj.alpha_fuel + obj.alpha_moderator)*(Tc - obj.Tc0) + obj.alpha_reflector*(Tr -obj.Tr0);
       
       end
       function dxdt = relative_neutron_flux(obj, t, x) 
           % DOCUMENT WHAT EACH VALUE REPRESENTS
           % x(1:10) = ni neutron flux of the nodes
           % x(11:70) = Concentrations of delayed groups for nodes
           % x(71:80) = ith Temperature of Fuel element nodes Tci
           % x(81:90) = ith Temperauter of Helium nodes Tdi
           Tr = obj.N*7+2*obj.N+1;% x(91) = Temperature of Reflector Tr
           Tu = obj.N*7+2*obj.N+2;% x(92) = Average Temperature in riser Tu
           Tlh = obj.N*7+2*obj.N+3;% x(93) = Temperature of lower helium header Tlh
           Toh = obj.N*7+2*obj.N+4;% x(94) = Temperature of outlet header Toh
           % x(95:104) = mass flowrate of helium at ith nodes h_mass
           Wu = obj.N*7+3*obj.N+5;% x(105) = mass flowrate of riser Wu
           Wlh = 145;% kg/s = mass flowrate fo lower helium header Wlh
           cores = obj.N*7+1; %This is the number to jump to core (fuel element) temperature
           downs = obj.N*7+obj.N+1; % This is the number to jump to downcomer (helium) temperature
           hmass = obj.N*7+2*obj.N+5; %number to get to masses  

           
           % Relative neutron flux for node 1
           rho_1 = obj.reactivity(x(cores),x(Tr)); % needs to be replaced with an equation
           
           dn1dt_term1 = (rho_1 - obj.beta - obj.lambda_ls_delayed_groups(1,1))/obj.lambda*x(1);
           dn1dt_term2 = (1/obj.lambda) * obj.coupling_coeffs_matrix(1,2) * x(2);
           dn1dt_term3 = obj.sum_beta_concentration_over_lambda(x(obj.N+1:obj.N+6), obj.lambda);
           dxdt(1) =  dn1dt_term1 + dn1dt_term2 + dn1dt_term3;
           
           % Relative neutron flux for ith nodes
           var = obj.N+7;
           for i = 2:(obj.N-1)
               
               rho_i = obj.reactivity(x(cores+i-1),x(Tr)); %need to be replaced withan equation
               
               dnidt_term1 = (rho_i - obj.beta - obj.coupling_coeffs_matrix(i,i))*(1/obj.lambda)*x(i);
               dnidt_term2 = (1/obj.lambda)*(obj.coupling_coeffs_matrix(i, i-1)*x(i-1) + obj.coupling_coeffs_matrix(i, i+1)*x(i+1));
               dnidt_term3 = obj.sum_beta_concentration_over_lambda(x(var:var+5), obj.lambda);
               dxdt(i) = dnidt_term1 + dnidt_term2 + dnidt_term3;
               var = var+6; 
           
           end
           % Relative neutron flux for node N
           rho_N = obj.reactivity(x(cores+obj.N-1),x(Tr)); % needs to be replaced with an equation
          
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
           %Values that need to be defined for now these are all random
           %number I created
           Kd = 116.9569; %W/(m^2*K) helium to the pellets/ pellets to helium 
           Ad = (4*3.14*(3.00E-02)^2)*420000/10; % m^2 volume of one sphere8 number of spheres/10 sections
           Kr = 80; % W/(m^2K) look into what the material is
           Ar = 207.35/10; % m^2 SA around column
           K = 100; % W/(m^2K)
           A = 28.27; % m^2 area of circle 
           Ku = 80; % W/(m^2K) look into what the material is, fluid and wall 
           Au = 241.9; %m^2 assuming reflector thickness = .5m 
           
           % Fuel Pile Temperature for 1st node
           dTc1dt_term1 = obj.P0*x(1);
           dTc1dt_term2 = - Kd*Ad*(x(cores)-x(downs));
           dTc1dt_term3 = - Kr*Ar*(x(cores)-x(Tr));
           dTc1dt_term4 = -K*A*(x(cores)-x(cores+1));
           dxdt(cores) = (dTc1dt_term1 +dTc1dt_term2 +dTc1dt_term3+ dTc1dt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);
           
           % Fuel Pile Temperature for ith node
           for i = 1:obj.N-2
               dTcidt_term1 = obj.P0*x(cores+i);
               dTcidt_term2 = - Kd*Ad*(x(cores+i)-x(downs+i));
               dTcidt_term3 = - Kr*Ar*(x(cores+i)-x(Tr));
               dTcidt_term4 = -K*A*(x(cores+i)-x(cores+i+1)) + K*A*(x(cores+i-1)- x(cores+i));
               dxdt(cores+i) = (dTcidt_term1 +dTcidt_term2 +dTcidt_term3 + dTcidt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);               
           end
           
           % Fuel Pile Temperature for Nth node
           dTcNdt_term1 = obj.P0*x(obj.N);
           dTcNdt_term2 = - Kd*Ad*(x(cores + obj.N-1)-x(downs + obj.N-1));
           dTcNdt_term3 = - Kr*Ar*(x(cores + obj.N-1)-x(Tr));
           dTcNdt_term4 = -K*A*(x(cores + obj.N-1 - 1)-x(cores+ + obj.N-1));
           dxdt(cores + obj.N-1) = (dTcNdt_term1 +dTcNdt_term2 +dTcNdt_term3+ dTcNdt_term4)/((1-obj.porosity)*obj.density_fuel*obj.volume_i*obj.C_fuel);           
           
           % Reflector Temperature           
           dTrdt_term1 = 0;
           for i = 0:obj.N -1
              dTrdt_term1 = dTrdt_term1 + Kr*Ar*(x(cores +i)-x(Tr)); 
           end
           dTrdt_term2 = -Ku*Au*(x(Tr)-x(Tu));
           dxdt(Tr) = (dTrdt_term1 + dTrdt_term2)/(obj.density_reflector*obj.volume_reflector*obj.C_reflector);
           
           % Downcomer Temperatrues
           %hmass node 1
           dxdt(hmass) = x(Wu) - x(hmass);
           %density node 2:N
           for i = 1:obj.N-1
               dxdt(hmass+i) = x(hmass+i-1) - x(hmass+i);
           end
           
           %Td node 1
           dTd1dt_term1 = Kd*Ad*(x(cores)- x(downs));
           dTd1dt_term2 = obj.Cp_helium*x(Wu)*(x(Tu)- x(downs));
           dxdt(downs) = (dTd1dt_term1 + dTd1dt_term2)/(obj.porosity*x(hmass)*obj.Cp_helium);
           
           %Td nodes 2:N
           for i = 1:obj.N-1
               dTdidt_term1 = Kd*Ad*(x(cores+i)-x(downs+i));
               dTdidt_term2 = obj.Cp_helium*x(hmass+i-1)*(x(downs+i-1)-x(downs+i));
               dxdt(downs+i)= (dTdidt_term1 + dTdidt_term2)/(obj.porosity*x(hmass+i)*obj.Cp_helium);
           end
           
           %Riser Temperature
           k = .01; %leakage ratio
           
           %hmass of riser Wu
           dxdt(Wu) = (1-k)*Wlh - x(Wu);
           
           %Tu 
           dTudt_term1 = obj.Cp_helium*(1-k)*Wlh*(x(Tlh) - x(Tu));
           dTudt_term2 = Ku*Au*(x(Tr) - x(Tu));
           dxdt(Tu) = (dTudt_term1 + dTudt_term2)/(x(Wu)*obj.Cp_helium);
           
           %Lower Header Tlh
           Tin = 250; % Helium Celsius 
           dxdt(Tlh) = (Tin - x(Tlh));
           
           %Outer Header Toh
           Woh = k*Wlh + x(hmass+obj.N-1);
           
           dTohdt_term1 = k*Wlh*(x(Tlh)-x(Toh));
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
            % This function is going to solve the ODE 
            tic
            tstep = .1;
            tspan_fix = tspan(1):tstep:tspan(2);
            [tout, x] = ode23tb(@obj.relative_neutron_flux, tspan_fix, x0);
            toc
       end

    end
    
end