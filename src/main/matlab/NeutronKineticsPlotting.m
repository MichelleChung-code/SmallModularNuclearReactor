% For plotting the neutron kinetics results 

% Plot neutron relative fluxxes
% TODO create labels dynamically, allow for different number of nodes

close all;

diff_plot_titles = ["Relative Neutron Fluxes of Nodes"];

for i = 1:N 
    diff_plot_titles(1 + i) = "Concentration Node " + num2str(i);
end

diff_plot_titles2 = ["Temperature of Fuel Elements", "Temperature of Helium", "Temperature of Reflector",...
"Average Temperature in Riser", "Temperature of Lower Helium Header",...
"Temperature of Outlet Header", "Nodal Mass Flow Rates of Helium",...
"Mass Flow Rate to Lower Helium Header", "Nodal Reactivities", "Rod Position - todo for future", "Nodal Power Output"];

diff_plot_titles = [diff_plot_titles, diff_plot_titles2];

diff_plot_ylabels = ["Neutron Fluxes",repelem(["Concentration"],...
    [N]),repelem(["Temperature"],[6]), repelem(["Mass Flow Rate"],[2]),"Reactivities", "Position", "Power Output (MW)"];


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
    ylabel(diff_plot_ylabels(i)), xlabel('Time, t')
    legend(legendInfo{starting_index:diff_plots_index_end(i)})
    starting_index = diff_plots_index_end(i) + 1;
    if i == 22 %todo dont hardcode
        total_at_end = num2str(round(sum(x(length(tout), 117:126))));
        annotation('textarrow',[0.5484 0.8778],[0.5775 0.5088],'String',strcat('Total Power:  ', total_at_end, 'MW'))
    end
end

disp("Plotting Completed");
