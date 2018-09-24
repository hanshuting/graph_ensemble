function gn_shuff_data(source_data, dest_folder, num_shuffles)
%GN_SHUFF_DATA Create shuffled copies of a dataset.
%
% Inputs
%   source_data: Filepath to .mat file containing dataset to shuffle in a variable named 'data'.
%   dest_folder: Path to folder shuffled copies will be created in.
%   num_shuffles (optional): Number of shuffled copies to create. Defaults to 100 if omitted.
if nargin < 3
    num_shuffles = 100;
end

addpath(genpath('..'));
load(source_data);
fprintf('Loaded: %s\n', source_data);
if exist('stimuli', 'var') ~= 1
    stimuli = [];
end
data_raw = data';
[~, filename, ~] = fileparts(source_data);
for i = 1:num_shuffles
        data = shuffle(data_raw,'exchange')';
        save(fullfile(dest_folder, ['shuffled_' filename '_' num2str(i) '.mat']),'data','stimuli');
end
fprintf('done shuffling data\n');