if(isempty(which('LoopyModelCollection')))
  addpath(genpath('~/src/crf/fwMatch'));
  exptPath = pwd;
  cd('~/src/crf/fwMatch');
  startup;
  cd(exptPath);  
end
dbclear all

%these are set in the script
if(isempty(who('SCRIPT')) || ...
   (~isempty(who('SCRIPT')) && ... 
    (isempty(who('LPATH')) || isempty(who('MODEL_TYPE')) || ...
     isempty(who('DEPTH')) || isempty(who('EXPT')))))
  LPATH = '~/src/crf/fwMatch/expt';
  MODEL_TYPE = '_tree';
  %ALL_EXPT = {'0', '45', '90', '135', '180', '225'};
  DEPTH = 1;
  EXPT = 'steph_0_01_tree';
  if(~isempty(who('SCRIPT')))
    warning(['You are calling from a script but not all vars are ' ...
             'set!!!']);
  end
end

SPATH = '~/data/steph_225_tree';
MODEL_PREFIX = 'steph_';
fload = sprintf('%s/%s/results/model_collection.mat', LPATH, ...
                EXPT);
load(fload);
fprintf('loaded model %s\n', fload);

best_model = model_collection.get_best_model;
inferred_model = best_model.inference_model;
clear model_collection best_model;
vars_N = size(inferred_model.x_train,2);


ss = sum(inferred_model.x_test,2);
ii = find(ss >= 8);
ii = reshape(ii, 1, length(ii));

for test_frame = ii
  infer_vars =  find(inferred_model.x_test(test_frame,:));
  observed_vars = setdiff(1:vars_N, infer_vars(1:end-2));
  [ll, qmap] = inferred_model.compute_conditional_likelihood_for_sample_true(inferred_model.x_test(test_frame,:), observed_vars);
  %disp(length(intersect(find(qmap), observed_vars))/length(observed_vars));
end

keyboard



CLL_0 = [];
CLL_1 = [];


for infer_id = 1:vars_N
  fprintf('inferring var %d\n', infer_id);
  observed_vars = setdiff(1:vars_N, infer_id);

  cond_l = inferred_model.compute_conditional_likelihood_all(observed_vars, ...
                                                    'train');
  CLL_0(infer_id,:) = cond_l.true_0;
  CLL_1(infer_id,:) = cond_l.true_1;
  %figure(1); plot(inferred_model.x_train(:,infer_id), 'go');
  %hold on;
  %plot(cond_l.true < -0.01, 'b+');
  %hold off;

end %all infer vars

save(['inferred_291114_luis_opto_ge3_' EXPT '_loopy.mat'], 'CLL_0', ...
     'CLL_1');

fprintf('saved: %s\n', ['inferred_291114_luis_opto_ge3_' EXPT '_loopy.mat']);




return

load ~/data/mouse/291114_luis_opto/ge3/ ...
     model_collection_001_s_tree.mat;

best_model = model_collection.get_best_model;
inferred_model = best_model.inference_model;
vars_N = size(inferred_model.x_train,2);

node_potentials = inferred_model.theta.node_potentials;
edge_potentials = inferred_model.theta.edge_potentials;
vars_N = size(inferred_model.x_train,2);
sample_count = 50000;
burn_in = 40000;
sample_interval = 200;
random_seed = 1;

[X, log_like] = gibbs_fast_sampler(node_potentials, edge_potentials, ...
                                   vars_N, sample_count, burn_in, ...
                                   sample_interval, random_seed);
marginals = sum(X)/sample_count;
plot(log_like);
X = logical(X);

log_like = log_like(vars_N*burn_in+1:end);

save('~/data/mouse/291114_luis_opto/ge3/001_s_tree_gibbs.mat', 'X', ...
     'sample_count', 'burn_in', 'sample_interval', 'random_seed', 'log_like');

%%%%
sampled_X = X(burn_in+1:end,:);
%ROIs:
load('~/data/mouse/291114_luis_opto/icaROIs_REG_C1V1GCaMP6s_m23_d1_001m_002m.mat');
load(['~/data/mouse/291114_luis_opto/' ...
      'SPK_REG_C1V1GCaMP6s_m23_d1_001m_002m.mat'], 'segcentroid');

label_positions(:,1) = segcentroid(:,1);
label_positions(:,2) = 240 - segcentroid(:,2);
for v = 1:vars_N
  variable_names{v} = sprintf('%d', v);
end
graph_to_dot(inferred_model.structure, 'filename', '001_s_tree.dot', ...
             'node_label', variable_names, 'positions', label_positions'.*2);

%neato -Tps -n  test_002.dot > test_002.ps; open test_002.ps



%%%%% 

fprintf('\nChow-liu test loglike: %f\n', ...
        chowliu_loglike(x_train,x_test));
fprintf('Bernoulli test loglike: %f\n', bernoulli_loglike(x_train, ...
                                                  x_test));


    model = treegmFit(x_train);
    avg_train_loglike = sum(treegmLogprob(model, x_train)) / ...
        size(x_train,1);
    avg_test_loglike = sum(treegmLogprob(model, x_test)) / ...
        size(x_test,1);
