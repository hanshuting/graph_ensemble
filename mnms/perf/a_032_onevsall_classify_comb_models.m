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
MODEL_TYPE = '_loopy';
MODEL_PREFIX = 'dff_mprs_';
ALL_EXPT = {'0', '45', '90', '135', '180', '225'}; %, '270', '315'};
DEPTH = 1;

ALL_LL = [];
for e = 1:length(ALL_EXPT)
  EXPT = sprintf('%s%s_%02d%s', MODEL_PREFIX, ALL_EXPT{e}, DEPTH, MODEL_TYPE);
  fload = sprintf('%s/LLs_%s_1vsall.mat',SPATH, EXPT);
  %save(fsave, 'LLs', 'EXPT');
  load(fload, 'LLs');

  fprintf('loaded model %s\n', fload);
  ALL_LL = [ALL_LL; LLs];
  clear   LLs;
end
LLs = ALL_LL;
clear ALL_LL;
models_N = length(ALL_EXPT);
fload = sprintf('%s/test_samples_%sd%02d.mat', SPATH, MODEL_PREFIX, ...
                DEPTH);
%save(fsave, 'test_samples', 'gnd');
load(fload, 'gnd');
samples_N = size(LLs,2);

%precision recall per class:
[vv, pred] = max(LLs, [], 1);

conf_matrix = zeros(models_N);
for m = 1:models_N
  for k = 1:models_N
    orient = str2num(ALL_EXPT{k});
    ii = intersect( find(pred == m), find(gnd == orient));
    conf_matrix(m,k) = length(ii)/length(find(gnd == orient));
  end
end
disp(conf_matrix);
conf_matrix_dup_orient = conf_matrix;

%0, 45, 90, 135, 180, 225
gnd((gnd == 180)) = 0;
gnd((gnd == 225)) = 45;
gnd((gnd == 270)) = 90;
gnd((gnd == 315)) = 135;
pred(pred == 5) = 1;
pred(pred == 6) = 2;
pred(pred == 7) = 3;
pred(pred == 8) = 4;
%this works because we just replaced the last two orientations:
models_N = length(unique(gnd));
conf_matrix = zeros(4);
for m = 1:models_N
  for k = 1:models_N
    orient = str2num(ALL_EXPT{k});
    ii = intersect( find(pred == m), find(gnd == orient));
    conf_matrix(m,k) = length(ii)/length(find(gnd == orient));
  end
end
disp(conf_matrix);


fsave = sprintf('%s/LLs_%sd%02d_1vsall%s.mat',SPATH, ...
                MODEL_PREFIX, DEPTH, MODEL_TYPE);
save(fsave, 'LLs', 'gnd', 'conf_matrix', 'ALL_EXPT', 'conf_matrix_dup_orient');
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