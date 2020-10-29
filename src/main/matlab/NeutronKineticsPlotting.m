% For plotting the neutron kinetics results 

% Plot neutron relative fluxxes
% TODO create labels dynamically, allow for different number of nodes

close all;

diff_plot_titles = ["Relative Neutron Fluxes of Nodes","Concentration Node 1",...
    "Concentration Node 2", "Concentration Node 3","Concentration Node 4",...
    "Concentration Node 5","Concentration Node 6","Concentration Node 7",...
    "Concentration Node 8","Concentration Node 9","Concentration Node 10",...
    "Temperature of Fuel Elements", "Temperature of Helium", "Temperature of Reflector",...
    "Average Temperature in Riser", "Temperature of Lower Helium Header",...
    "Temperature of Outlet Header", "Nodal Mass Flow Rates of Helium",...
    "Mass Flow Rate to Lower Helium Header"];

diff_plot_ylabels = ["Neutron Fluxes",repelem(["Concentration"],...
    [10]),repelem(["Temperature"],[6]), repelem(["Mass Flow Rate"],[2])];
diff_plots_index_end = [10, 11+5*(1), 10+6*(2), 10+6*(3), 10+6*(4),...
    10+6*(5), 10+6*(6), 10+6*(7),10+6*(8), 10+6*(9), 10+6*(10), 80, 90, 91,...
    92 ,93 ,94,  104, 105];
series_num = length(x);
legendInfo = cell(1, series_num); % predefine size for performance

const_label = "Node ";
delayed_groups = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
for i = 1:series_num 
    if (1 <= i) && (i <= 10); legendInfo{i} = [const_label + num2str(i)]; end
    if (10 <= i) && (i <= 70)   
        n = i-11;
        if ismember(11+ n*6, delayed_groups*6+11); const_label = "Delayed Group "; end
        group_num = i-([delayed_groups]*6+11);
        legendInfo{i} = [const_label + num2str(min(group_num(group_num>=0))+1)];
    end
    if (71 <= i) && (i <= 80); legendInfo{i} = ["Node " + num2str(i-70)]; end
    if (81 <= i) && (i <= 90); legendInfo{i} = ["Node " + num2str(i-80)]; end
    if (91 <= i) && (i <= 94); legendInfo{i} = ["Temperature"]; end
    if (95 <= i) && (i <= 104); legendInfo{i} = ["Node "+ num2str(i-94)]; end
    if (i == 105); legendInfo{i} = ["Mass Flow Rate"]; end
end

starting_index = 1;
for i=1:length(diff_plot_titles)
    figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
    title(diff_plot_titles(i))
    ylabel(diff_plot_ylabels(i)), xlabel('Time, t')
    legend(legendInfo{starting_index:diff_plots_index_end(i)})
    starting_index = diff_plots_index_end(i) + 1;
end

