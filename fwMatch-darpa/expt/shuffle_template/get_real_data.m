function [data,variable_names] = get_real_data()

%load('/vega/brain/users/sh3276/data/luis/m24_d1_05_opto_on.mat');
load(['/vega/brain/users/sh3276/data/luis/m24_vis.mat']);
fprintf('Loaded: %s\n', ['/vega/brain/users/sh3276/data/luis/m24_vis.mat']);

%data is time_frames by number_of_neurons; it is logical andsparse
%Kate
%data = ALL_DATA' > 0.3;
data = full(data);
N = size(data,2);
variable_names = {};
for i = 1:N
    variable_names(end+1) = {int2str(i)};
end
fprintf('data is : %d, %d\n', size(data,1), size(data,2));
end
