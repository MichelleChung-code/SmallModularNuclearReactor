% Required Data
% - D is the diffusion coefficient - N
% - V is the nodal volume - Y
% - d is the nodal height - Y 
% - sigma_a per ith node = macroabsorption cross section of node i - N
% - sigma_f per fth node = fission cross section of node i - N
% - v is the fission number - N 

% Calculate diffusion coefficient ourselves from SS values and assume it
% remains constant
% Calculated following: 
% https://www.nuclear-power.net/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/
clear, clc
uranium_mass_number = 235;
mu_avg = 2/(3*uranium_mass_number);

% scattering cross section data from: https://www.osti.gov/etdeweb/servlets/purl/20956173
% file located in "C:\Users\tkdmc\OneDrive - University of Calgary\Capstone_Group25_CHEMENGG\Reactor_Modelling\Reactor_Core_Modelling\Michelle_Research_Winter_2021\Cross_Sections_And_Neutron_Yields_for_U235.pdf"
sigma_s = 15E-24; % 10 +- 2 barns, 1 barn = 10^-24 cm^2 scattering cross section or microscopic cross-section 
uranium_density = 19.1; % g/cm^3
avogadro_num = 6.0221409E23; % nuclei/mol
atomic_number_density = uranium_density * avogadro_num / uranium_mass_number;  % nuclei / cm^3
macroscopic_cross_section = sigma_s * atomic_number_density; % cm^-1
D = 1/(3*macroscopic_cross_section*(1-mu_avg)); % Diffusion Coefficient [cm]
transport_mean_free_path = 3*D; 

% literature sources from the tables were arbitrarily chosen 
absorption_cross_section = 695E-24; % cm^2 take value from Table 3, Nikitan et al. (1956)
fission_cross_section = 605E-24; % cm^2 take value from Table 4, Saplakoglu (1958) 
fission_number = 2.5; % https://www.world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/physics-of-nuclear-energy.aspx#:~:text=Fission%20of%20U%2D235%20nuclei,absorbed%20in%20non%2Dfission%20reactions. hopefully they are referring to the number of neutrons released per fission...

% Nodal Items 
reactor_core_diameter = 3; % m
reactor_core_height = 11; % m 

num_nodes = 10;
node_height = reactor_core_height/num_nodes; % m
node_diameter = reactor_core_diameter; % m 

A = pi*(node_diameter/2)^2; %m^2
V = A*node_height; %m^3


a_matrix = zeros(num_nodes, num_nodes);

% Everything in SI units
D = D*0.01; % m diffusion coefficient
absorption_cross_section = absorption_cross_section * 0.0001; %m^2
fission_cross_section = fission_cross_section* 0.0001; %m^2

% for initial neutron fluxes (use SS values, which are the initial
% conditions)
% Therefore, the maximum amount of nodes we can have is 10 

from_csv = readtable('Initial Values.csv');
nodal_neutron_fluxes_0 = table2array(from_csv(1:10,4));

for i=1:num_nodes
    if i==1
        j = [2];
        a_matrix(i,i) = D*A/(fission_number*fission_cross_section*V*node_height);
    elseif i == num_nodes
        j = [num_nodes-1];
        a_matrix(i,i) = (D*A)/(fission_number*fission_cross_section*V*node_height);
    else
        j = [i-1, i+1];
        a_matrix(i,i) = (D*A*(2*(1/fission_cross_section)))/(fission_number*V*node_height);
    end
    
    for k=1:length(j)
        a_matrix(i,j(k)) = (nodal_neutron_fluxes_0(j(k))/nodal_neutron_fluxes_0(i)) * (D*A/(fission_number*fission_cross_section*V*node_height));
    end
    
end

common_factor = 1E-25;

a_matrix = a_matrix.*common_factor; % okay because we normalize to node 1 prior to plotting anyways.  This just makes the numbers easier to deal with.

disp(a_matrix)
