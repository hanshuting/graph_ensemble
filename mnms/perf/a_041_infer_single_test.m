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
  DATASET = 'test';

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

CLL_0 = []; CLL_1 = [];
test_frames_N = size(inferred_model.x_test,1);

fload = sprintf('%s/INFERRED/%s_%s.mat', SPATH, EXPT, DATASET);
%save(fsave, 'CLL_0', 'CLL_1');
if(exist(fload, 'file'))
  load(fload);
  fprintf('LOADED: %s\n', fload);
  last_frame = size(CLL_0,2);
  frames_run = (last_frame-1):test_frames_N;
  fprintf('running from frame: %d\n', last_frame-1);
else
  frames_run = 1:test_frames_N;
end

fsave = sprintf('%s/INFERRED/%s_%s.mat', SPATH, EXPT, DATASET);
    
for frame = frames_run
  sample_orig = inferred_model.x_test(frame,:);
  if(mod(frame, 50) == 0)
    fprintf('%d/%d\n', frame, test_frames_N);
    save(fsave, 'CLL_0', 'CLL_1');
    fprintf('saved: %s\n', fsave);
    disp(datestr(now));
  end
  for infer_id = 1:vars_N
    %fprintf('inferring var %d\n', infer_id);
    observed_vars = setdiff(1:vars_N, infer_id);
    sample = sample_orig;
    sample(infer_id) = 0;
    if(strcmp(MODEL_TYPE, 'loopy') == 1)
      cond_l = inferred_model.compute_conditional_likelihood_for_sample_approx(sample, ...
                                                        observed_vars);
      CLL_0(infer_id,frame) = cond_l;
      sample(infer_id) = 1;      
      cond_l = inferred_model.compute_conditional_likelihood_for_sample_approx(sample, ...
                                                        observed_vars);
      CLL_1(infer_id,frame) = cond_l;

    else
      cond_l = inferred_model.compute_conditional_likelihood_for_sample_true(sample, ...
                                                        observed_vars);
      CLL_0(infer_id,frame) = cond_l;    
      sample(infer_id) = 1;
      cond_l = inferred_model.compute_conditional_likelihood_for_sample_true(sample, ...
                                                        observed_vars);
      CLL_1(infer_id,frame) = cond_l;    
    end
  end
end %all infer vars

fsave = sprintf('%s/INFERRED/%s_%s.mat', SPATH, EXPT, DATASET);
save(fsave, 'CLL_0', 'CLL_1');
fprintf('saved: %s\n', fsave);
disp(datestr(now));
return

