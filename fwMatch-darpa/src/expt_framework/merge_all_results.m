function [model_collection] = merge_all_results(results_folder_path)
% MERGE_ALL_RESULTS
%
% Input
%   results_folder_path (Optional): String to folder containing result*.mat
%       files to merge. MUST include a trailing /.
    if nargin < 1
        results_folder_path = pwd;
    end

    % iterate through all result files merging the model collections
    result_files = dir(sprintf('%sresult*.mat', results_folder_path));
    for i = 1:length(result_files)
        split_result = load(sprintf('%s%s', results_folder_path, result_files(i).name), 'model_collection');
        if i == 1
            merged_model_collection = split_result.model_collection;
        else
            merged_model_collection = merged_model_collection.merge_model_collections(split_result.model_collection);
        end
    end
    model_collection = merged_model_collection;

    save(sprintf('%smodel_collection.mat', results_folder_path), 'model_collection', '-v7.3');
end
