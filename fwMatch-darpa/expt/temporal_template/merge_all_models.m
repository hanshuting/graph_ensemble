
exptPath = pwd;
cd('../../');
% startup;

% iterate through all result files merging the model collections
result_files = dir(sprintf('%s/results/result*.mat', exptPath));
for i = 1:length(result_files)
    split_result = load(sprintf('%s/results/%s', exptPath, result_files(i).name), 'model_collection');
    if i == 1
        merged_model_collection = split_result.model_collection;
    else
        merged_model_collection = merged_model_collection.merge_model_collections(split_result.model_collection);
    end   
end
model_collection = merged_model_collection;

save(sprintf('%s/results/model_collection.mat', exptPath), 'model_collection', '-v7.3');
%fsave = exptPath(strfind(exptPath, 'expt/')+length('expt/'):end);
%save(sprintf('~/data/steph_225_tree/model_collection_%s.mat', fsave), 'model_collection', '-v7.3');

cd(exptPath);
