
[~,savename,~] = fileparts(pwd);
savepath = '/vega/brain/users/sh3276/results/luis/models/';

load('./results/model_collection.mat');
best_model = getBestModel(model_collection);

save([savepath savename '_best_model_org.mat'],'best_model','-v7.3');
