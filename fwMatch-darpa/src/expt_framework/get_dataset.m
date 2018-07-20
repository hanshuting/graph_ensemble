function [data,variable_names, stimuli] = get_dataset(data_filepath)
load(data_filepath);
fprintf('Loaded: %s\n', data_filepath);
%data is time_frames by number_of_neurons
data = full(data);
N = size(data,2);
fprintf('data is : %d, %d\n', size(data,1), size(data,2));
variable_names = {};
for i = 1:N
    variable_names(end+1) = {int2str(i)};
end
if exist('stimuli', 'var') == 1
    stimuli = full(stimuli);
    fprintf('stimuli is : %d, %d\n', size(stimuli,1), size(stimuli,2));
else
    stimuli = [];
end
end
