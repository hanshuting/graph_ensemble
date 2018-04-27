function [data,variable_names] = get_real_data()
load(['/Users/jonathanshor/GitHub/graph_ensemble/data' ...
          '/shuffled_m52_d1_opto_high_all_add_neuron_7.mat']);
fprintf('Loaded: %s\n', ['/Users/jonathanshor/GitHub/graph_ensemble/data' ...
'/shuffled_m52_d1_opto_high_all_add_neuron_7.mat']);
%data is time_frames by number_of_neurons
data = full(data);
N = size(data,2);
variable_names = {};
for i = 1:N
	variable_names(end+1) = {int2str(i)};
end
fprintf('data is : %d, %d\n', size(data,1), size(data,2));
end
