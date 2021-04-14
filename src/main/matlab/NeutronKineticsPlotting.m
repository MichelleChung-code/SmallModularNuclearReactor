% For plotting the neutron kinetics results 

% Plot neutron relative fluxxes

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
global reactivity_control_rod_results
global time_plot

starting_index = 1;
for i=1:length(diff_plot_titles)
%     hard coded y-limits for producing steady state plots
%     if i == 14 + 1
%         xlim([100 inf]) 
%         ylim([250 255])
%     elseif i == 2 + 1
%         xlim([100 inf]) 
%         ylim([6.9e-4 7e-4])
%     elseif i == 4 + 1
%         xlim([100 inf]) 
%         ylim([3.3e-3 3.4e-3])
%     elseif i == 3 + 1
%         xlim([100 inf]) 
%         ylim([1.2e-3 1.3e-3])
%     elseif i == 15 + 1 
%         xlim([100 inf]) 
%         ylim([250 255])
%     elseif i == 17 + 1 
%         xlim([100 inf]) 
%         ylim([740 770])
%     elseif i == 5 + 1 
%         xlim([100 inf]) 
%         ylim([0.01 0.02])
%     elseif i == 6 + 1 
%         xlim([100 inf]) 
%         ylim([0.04 0.05])  
%     elseif i == 7 + 1 
%         xlim([100 inf]) 
%         ylim([0.18 0.19])
%     elseif i == 8 + 1 
%         xlim([100 inf]) 
%         ylim([0.8 0.9])
%     elseif i == 9 + 1 
%         xlim([100 inf]) 
%         ylim([1.8 1.9])
%     elseif i == 10 + 1 
%         xlim([100 inf]) 
%         ylim([1.2 1.3])
%     elseif i == 11 + 1 
%         xlim([100 inf]) 
%         ylim([0.6 0.7])
%     end
%     
    
    
    
    
    if diff_plot_titles(i) == "Nodal Reactivity"
%         figure(i);
%         hold on;
%         plot(time_plot, reactivity_control_rod_results(:, 1:N)), grid on
%         for j=1:N
%             coefficients = polyfit(time_plot, reactivity_control_rod_results(:, j), 6);
%             xFit = linspace(min(time_plot), max(time_plot), 1000);
%             yFit = polyval(coefficients , xFit);
%             plot(xFit, yFit, 'r-'); % Plot fitted line.
%         end
%         hold off;

% start x axis at 10 seconds for the derivative approx
        figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
        xlim([100 inf])         
        
    elseif diff_plot_titles(i) == "Rod Position"
%         figure(i);
%         hold on;
%         plot(time_plot, reactivity_control_rod_results(:, N + 1)), grid on    
%         coefficients = polyfit(time_plot, reactivity_control_rod_results(:, N + 1), 5);
%         xFit = linspace(min(time_plot), max(time_plot), 1000);
%         yFit = polyval(coefficients , xFit);
%         plot(xFit, yFit, 'r-'); % Plot fitted line.
%         hold off;
% start x axis at 10 seconds for the derivative approx
        figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
        xlim([1000 inf]) 
        ylim([0 10])
    else
        figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
    end
    
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

figure(i+2), plot(1:N, reactivity_control_rod_results(end,1:N));
title("Nodal Final Reactivity");
ylabel("Reactivity (MW)"), xlabel("Node Number");

figure(i+3), plot(1:N, x(end,1:N));
title("Nodal Final Relative Neutron Flux");
ylabel("Relative Neutron Flux"), xlabel("Node Number");
disp("Plotting Completed");
