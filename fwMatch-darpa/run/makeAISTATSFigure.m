%% Generating numbers and a pretty figure for the AISTATS 2016 submission.

% First, load the data and get the best model.
%load expt/loopy_model_collection_real_data/results_2012_2014/model_collection.mat
load expt/loopy_model_collection_real_data/results_10000/model_collection.mat
%load expt/loopy_model_collection_real_data/results_20000/model_collection.mat
model_collection.get_best_model()

%% Compute baseline likelihoods
chowliu_loglike(model_collection.x_train,model_collection.x_test)
bernoulli_loglike(model_collection.x_train,model_collection.x_test)

%% Get the labels and the figures
load expt/loopy_model_collection_real_data/variables_short.mat

%full_variables_st = load('expt/loopy_model_collection_real_data/variables_original.mat');
%full_variables = full_variables_st.variable_names;

% Filter the column names
variable_names = cellfun(@(x) strrep(x, 'diff_', 'd'), variable_names, 'UniformOutput', false);
variable_names = cellfun(@(x) strrep(x, '%_', '%d'),   variable_names, 'UniformOutput', false);
variable_names = cellfun(@(x) strrep(x, '_WORLD', ''), variable_names, 'UniformOutput', false);

%% Finally, visualize the graph
% Assign
model_collection.variable_names = variable_names;
model_collection.visualize_best_model_graph()
