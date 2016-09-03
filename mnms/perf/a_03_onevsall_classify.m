if(isempty(which('LoopyModelCollection')))
  addpath(genpath('~/src/crf/fwMatch'));
  exptPath = pwd;
  cd('~/src/crf/fwMatch');
  startup;
  cd(exptPath);  
end
dbclear all


LPATH = '~/src/crf/fwMatch/expt';
SPATH = '~/data/steph_225_tree';
MODEL_TYPE = '_tree';
MODEL_PREFIX = 'steph_';
ALL_EXPT = {'0', '45', '90', '135', '180', '225'};

DEPTH = 1;

for e = 1:length(ALL_EXPT)
  EXPT = ALL_EXPT{e};
  fload = sprintf('%s/%s%s_%02d%s/results/model_collection.mat', LPATH, ...
                  MODEL_PREFIX, EXPT, DEPTH, MODEL_TYPE);
  %load(fload);
  fload = sprintf('%s/LLs_%s_1vsall.mat',SPATH, EXPT);
  %save(fsave, 'LLs', 'EXPT');

  fprintf('loaded model %s\n', fload);

  best_model = model_collection.get_best_model;
  inferred_model = best_model.inference_model;
  clear best_model;
  vars_N = size(inferred_model.x_train,2);

  %try to reduce mem usage:
  inferred_model.x_train = [];

  MODELS(e) = inferred_model;
  clear inferred_model model_collection;
end


models_N = length(MODELS);
LLs = zeros(models_N);


fload = sprintf('%s/test_samples_d%02d.mat', SPATH, DEPTH);
%save(fsave, 'test_samples', 'gnd');
load(fload);

%test_samples = [];
%gnd = [];
%for test_model = 1:models_N
%  test_samples = [test_samples; MODELS(test_model).x_test];
%  gnd = [gnd, test_model.*ones(1,size(MODELS(test_model).x_test,1))];
%end

samples_N = size(test_samples,1);

%keyboard

for train_model = 1:models_N
  fprintf('computing for train model %d\n', train_model);
  %for test_model = 1:models_N
  for s = 1:samples_N
    LLs(train_model, s) = compute_avg_log_likelihood(MODELS(train_model).theta.node_potentials, ...
                                                     MODELS(train_model).theta.edge_potentials, ...
                                                     MODELS(train_model).theta.true_logZ, test_samples(s,:));
  end
end




%precision recall per class:

[vv, pred] = max(LLs, [], 1);

conf_matrix = zeros(models_N);
for m = 1:models_N
  for k = 1:models_N
    ii = intersect( find(pred == m), find(gnd == k));
    conf_matrix(m,k) = length(ii)/length(find(gnd == k));
  end
end


disp(conf_matrix);

fsave = sprintf('%s/LLs_%sonevsall%s.mat',SPATH, MODEL_PREFIX, MODEL_TYPE);
save(fsave, 'LLs', 'gnd', 'conf_matrix', 'ALL_EXPT');
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