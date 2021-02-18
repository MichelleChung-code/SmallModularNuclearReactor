% For plotting the neutron kinetics results 

% Plot neutron relative fluxxes
% TODO create labels dynamically, allow for different number of nodes

close all;
deg_sign = char(0176);
diff_plot_titles = ["Relative Nodal Neutron Flux"];

for i = 1:N 
    diff_plot_titles(1 + i) = "Relative Delayed Neutron Concentration of Node " + num2str(i);
end

diff_plot_titles2 = ["Temperature of Fuel Elements", "Temperature of Helium in Nodes", "Temperature of Reflector",...
"Average Temperature in Riser", "Temperature of Lower Helium Header",...
"Temperature of Outlet Header", "Mass Flow Rates of Helium through Nodes",...
"Mass Flow Rate to Lower Helium Header", "Nodal Reactivity", "Rod Position", "Nodal Power Output"];

diff_plot_titles = [diff_plot_titles, diff_plot_titles2];

diff_plot_ylabels = ["Relative Neutron Flux",repelem(["Relative Concentration"],...
    [N]),repelem([strcat("Temperature (", deg_sign, "C)")],[6]), repelem(["Mass Flow Rate (kg/s)"],[2]),"Reactivity", "Position", "Power Output (MW)"];


for i = 1:N+1
   diff_plots_index_end(i) = N + 6*(i-1);
end

diff_plots_index_end2 = [N*8, 9*N, 9*N + 1,...
    9*N + 2 ,9*N + 3 ,9*N + 4,  10*N + 4, 10*N + 5, 11*N + 5, 11*N + 6,12*N + 6];

diff_plots_index_end = [diff_plots_index_end, diff_plots_index_end2];

series_num = length(x);
legendInfo = cell(1, series_num); % predefine size for performance

const_label = "Node ";
delayed_groups = 0:1:N;
for i = 1:series_num 
    if (1 <= i) && (i <= N); legendInfo{i} = [const_label + num2str(i)]; end
    if (N < i) && (i <= 7*N)   
        n = i-(N+1);
        if ismember((N+1)+ n*6, delayed_groups*6+(N+1)); const_label = "Delayed Group "; end
        group_num = i-([delayed_groups]*6+(N+1));
        legendInfo{i} = [const_label + num2str(min(group_num(group_num>=0))+1)];
    end
    if (7*N + 1 <= i) && (i <= 8*N); legendInfo{i} = ["Node " + num2str(i-(7*N))]; end
    if (8*N + 1 <= i) && (i <= 9*N); legendInfo{i} = ["Node " + num2str(i-(8*N))]; end
    if (9*N + 1 <= i) && (i <= 9*N + 4); legendInfo{i} = ["Temperature"]; end
    if (9*N + 5 <= i) && (i <= 10*N + 4); legendInfo{i} = ["Node "+ num2str(i-(9*N + 4))]; end
    if (i == 10*N + 5); legendInfo{i} = ["Mass Flow Rate"]; end
    if (10*N + 5 < i) && (i <= 11*N + 5); legendInfo{i} = ["Node "+ num2str(i-(10*N + 5))]; end
    if (i == 11*N + 6); legendInfo{i} = ["Rod Position"]; end
    if (11*N + 6 < i) && (12*N + 6); legendInfo{i} = ["Node "+ num2str(i-(11*N + 6))]; end
end

starting_index = 1;
for i=1:length(diff_plot_titles)
    figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
    title(diff_plot_titles(i))
    ylabel(diff_plot_ylabels(i)), xlabel('Time (s)')
    legend(legendInfo{starting_index:diff_plots_index_end(i)}) %'Location','southeast'
    starting_index = diff_plots_index_end(i) + 1;
    
    [num_row_x, num_col_x] = size(x);
    total_graphs = 22 - (10-N); % 22 is the number of graphs for 10 nodes 
    if i == total_graphs 
        total_at_end = num2str(round(sum(x(length(tout), num_col_x - (N-1): num_col_x))));
        annotation('textbox',[.4 .5 .3 .4],'String',strcat('Total Power:  ', total_at_end, 'MW'),'FitBoxToText','on')
    end
end
global num_col_x_without_PI

figure(i+1), plot(1:N, x(end,num_col_x_without_PI + 1:num_col_x_without_PI + N));
title("Nodal Final Power Output");
ylabel("Power Output (MW)"), xlabel("Node Number");

figure(i+2), plot(1:N, x(end,num_col_x_without_PI - N :num_col_x_without_PI - 1));
title("Nodal Final Reactivity");
ylabel("Reactivity (MW)"), xlabel("Node Number");

figure(i+3), plot(1:N, x(end,1:N));
title("Nodal Final Relative Neutron Flux");
ylabel("Relative Neutron Flux"), xlabel("Node Number");
disp("Plotting Completed");
