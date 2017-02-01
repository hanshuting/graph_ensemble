function [data,variable_names] = get_real_data()
load(['/vega/brain/users/sh3276/data/luis' ...
          '/m21_d2_vis_all_high_add_neuron.mat']);
fprintf('Loaded: %s\n', ['/vega/brain/users/sh3276/data/luis' ...
'/m21_d2_vis_all_high_add_neuron.mat']);
%data is time_frames by number_of_neurons
data = full(data);
N = size(data,2);
variable_names = {};
for i = 1:N
	variable_names(end+1) = {int2str(i)};
end
fprintf('data is : %d, %d\n', size(data,1), size(data,2));
end
