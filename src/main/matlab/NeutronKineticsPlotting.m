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
"Mass Flow Rate to Lower Helium Header", "Nodal Reactivities"];

diff_plot_titles = [diff_plot_titles, diff_plot_titles2];

diff_plot_ylabels = ["Neutron Fluxes",repelem(["Concentration"],...
    [N]),repelem(["Temperature"],[6]), repelem(["Mass Flow Rate"],[2]),"Reactivities"];


for i = 1:N+1
   diff_plots_index_end(i) = N + 6*(i-1);
end

diff_plots_index_end2 = [N + N*6 + N, N + N*6 + 2*N, N + N*6 + 2*N + 1,...
    N + N*6 + 2*N + 2 ,N + N*6 + 2*N + 3 ,N + N*6 + 2*N + 4,  N + N*6 + 3*N + 4, N + N*6 + 3*N + 5, 11*N + 5];

diff_plots_index_end = [diff_plots_index_end, diff_plots_index_end2];

series_num = length(x);
legendInfo = cell(1, series_num); % predefine size for performance

const_label = "Node ";
delayed_groups = 0:1:N;
for i = 1:series_num 
    if (1 <= i) && (i <= N); legendInfo{i} = [const_label + num2str(i)]; end
    if (N < i) && (i <= N+6*N)   
        n = i-(N+1);
        if ismember((N+1)+ n*6, delayed_groups*6+(N+1)); const_label = "Delayed Group "; end
        group_num = i-([delayed_groups]*6+(N+1));
        legendInfo{i} = [const_label + num2str(min(group_num(group_num>=0))+1)];
    end
    if (N+6*N + 1 <= i) && (i <= N + N*6 + N); legendInfo{i} = ["Node " + num2str(i-(N+6*N))]; end
    if (N + N*6 + N + 1 <= i) && (i <= N + N*6 + 2*N); legendInfo{i} = ["Node " + num2str(i-(N + N*6 + N))]; end
    if (N + N*6 + 2*N + 1 <= i) && (i <= N + N*6 + 2*N + 4); legendInfo{i} = ["Temperature"]; end
    if (N + N*6 + 2*N + 5 <= i) && (i <= N + N*6 + 3*N + 4); legendInfo{i} = ["Node "+ num2str(i-(N + N*6 + 2*N + 4))]; end
    if (i == N + N*6 + 3*N + 5); legendInfo{i} = ["Mass Flow Rate"]; end
    if (N + N*6 + 3*N + 5 < i) && (i <= 11*N + 5); legendInfo{i} = ["Node "+ num2str(i-(10*N + 5))]; end
end

starting_index = 1;
for i=1:length(diff_plot_titles)
    figure(i), plot(tout, x(:, starting_index:diff_plots_index_end(i))), grid on
    title(diff_plot_titles(i))
    ylabel(diff_plot_ylabels(i)), xlabel('Time, t')
    legend(legendInfo{starting_index:diff_plots_index_end(i)})
    starting_index = diff_plots_index_end(i) + 1;
end

