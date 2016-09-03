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
  MODEL_PREFIX = 'rawb_';
  if(~isempty(who('SCRIPT')))
    warning(['You are calling from a script but not all vars are ' ...
             'set!!!']);
  end
end
SPATH = '~/data/steph_225_tree';
fload = sprintf('%s/%s/results/model_collection.mat', LPATH, ...
                EXPT);
load(fload);
fprintf('loaded model %s\n', fload);

best_model = model_collection.get_best_model;
inferred_model = best_model.inference_model;
clear model_collection best_model;
vars_N = size(inferred_model.x_train,2);

node_pot = inferred_model.theta.node_potentials;
edge_pot = inferred_model.theta.edge_potentials;
if(isfield(inferred_model.theta, 'true_logZ') )
  logZ = inferred_model.theta.true_logZ;
else
  logZ = inferred_model.theta.logZ;
end
clear inferred_model;


%load test data for this depth:
fload = sprintf('%s/test_samples_%sd%02d.mat', SPATH, MODEL_PREFIX, DEPTH);
%save(fsave, 'test_samples', 'gnd');
load(fload);
fprintf('LOADED: %s\n', fload);
LLs = zeros(1,length(gnd));
samples_N = size(test_samples,1);

if(islogical(test_samples) || length(unique(test_samples(:))) == 2)
  binary_model = 1;
else
  binary_model = 0;
end

for s = 1:samples_N
  if(binary_model == 1)
    LLs(s) = compute_avg_log_likelihood(node_pot, edge_pot, logZ, ...
                                        test_samples(s,:));
  else
    LLs(s) = compute_avg_log_likelihood_hidden(node_pot, edge_pot, logZ, ...
                                        test_samples(s,:));
  end
  %inferred_model.theta.node_potentials, ...
  %inferred_model.theta.edge_potentials, ...
  %inferred_model.theta.true_logZ, test_samples(s,:));
  if(mod(s,500) == 0)
    fprintf('done test frame %d/%d\n', s, samples_N);
  end
end
fsave = sprintf('%s/LLs_%s_1vsall.mat',SPATH, EXPT);
save(fsave, 'LLs', 'EXPT');
fprintf('SAVED: %s\n', fsave);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%classification per class:

[vv, pred] = max(LLs, [], 1);

conf_matrix = zeros(models_N);
for m = 1:models_N
  for k = 1:models_N
    ii = intersect( find(pred == m), find(gnd == k));
    conf_matrix(m,k) = length(ii)/length(find(gnd == k));
  end
end


disp(conf_matrix);

fsave = sprintf('%s/LLs_%s_1vsall.mat',SPATH, EXPT);
save(fsave, 'LLs', 'EXPT');
fprintf('SAVED: %s\n', fsave);


return

TRAIN_TPR = sparse([]); %size(inferred_model.x_train));
TRAIN_TNR = sparse([]); %size(inferred_model.x_train));
TRAIN_ACC = sparse([]); %size(inferred_model.x_train));

TEST_TPR = sparse([]); %size(inferred_model.x_test));
TEST_TNR = sparse([]); %size(inferred_model.x_test));
TEST_ACC = sparse([]); %size(inferred_model.x_test));


%infer_id = 40;

fprintf('debug\n');
keyboard

for infer_id = 40
  %1:vars_N
  fprintf('\ninferring var %d\n', infer_id);
  observed_vars = setdiff(1:vars_N, infer_id);

  cond_l = inferred_model.compute_conditional_likelihood_all(observed_vars, 'train');
  %figure(1); plot(inferred_model.x_train(:,infer_id), 'go');
  %hold on;
  %plot(cond_l.true < -0.01, 'b+');
  %hold off;

  gnd = (inferred_model.x_train(:,infer_id) == 1)';
  pred = cond_l.true < -1;

  %true positive:
  gnd_pos = find(gnd == 1);
  pred_pos = find(pred == 1);
  tp = length(intersect(gnd_pos, pred_pos)); 
  %true negative:
  gnd_neg = find(gnd == 0);
  pred_neg = find(pred == 0);
  tn = length(intersect(gnd_neg, pred_neg)); 
  %accuracy: (tp + tn)/(p+n) = (tp+tn)/N
  accuracy = (tp+tn)/length(gnd);

  fprintf('TRAIN: tp: %0.2f, tn: %0.2f, acc: %0.2f\n', ...
          tp, tn, accuracy);

  TRAIN_TPR(infer_id) = tp/length(gnd_pos);
  TRAIN_TNR(infer_id) = tn/length(gnd_neg);
  TRAIN_ACC(infer_id) = accuracy;


  cond_l = inferred_model.compute_conditional_likelihood_all(observed_vars, 'test');
  %figure(2); plot(inferred_model.x_test(:,infer_id), 'go');
  %hold on;
  %plot(cond_l.true < -0.01, 'b+');
  %hold off;

  gnd = (inferred_model.x_test(:,infer_id) == 1)';
  pred = cond_l.true < -1;

  %true positive:
  gnd_pos = find(gnd == 1);
  pred_pos = find(pred == 1);
  tp = length(intersect(gnd_pos, pred_pos)); 
  %true negative:
  gnd_neg = find(gnd == 0);
  pred_neg = find(pred == 0);
  tn = length(intersect(gnd_neg, pred_neg)); 
  %accuracy: (tp + tn)/(p+n) = (tp+tn)/N
  accuracy = (tp+tn)/length(gnd);

  fprintf('TEST: tp: %0.2f, tn: %0.2f, acc: %0.2f\n', ...
          tp, tn, accuracy);
  TEST_TPR(infer_id) = tp/length(gnd_pos);
  TEST_TNR(infer_id) = tn/length(gnd_neg);
  TEST_ACC(infer_id) = accuracy;

end %all infer vars


fsave = ['~/data/mouse/291114_luis_opto/ge3/inferred_results_' EXPT ...
      '.mat'];

save(fsave, 'TRAIN_TPR', 'TRAIN_TNR', 'TRAIN_ACC', ...
     'TEST_TPR', 'TEST_TNR', 'TEST_ACC');
fprintf('SAVED: %s\n', fsave);