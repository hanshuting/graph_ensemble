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
MODEL_PREFIX = 'steph_';


ORIENTS = (0,45,90,135,180,225,270,315);

count = 1;
ACC = [];
for ORIENT = ORIENTS
  EXPT = sprintf('%s%d_%02d', MODEL_PREFIX,ORIENT, DEPTH, MODEL_TYPE);
  fload = sprintf('%s/INFERRED/%s_%s.mat', EXPT, DATASET);
  %these are vars by frames:
  %save(fsave, 'CLL_0', 'CLL_1');
  load(fload)
  fprintf('loaded: %s\n', fload);

  fload = sprintf('%s/%s/results/model_collection.mat', LPATH, ...
                  EXPT);
  load(fload);
  fprintf('loaded model %s\n', fload);

  best_model = model_collection.get_best_model;
  inferred_model = best_model.inference_model;
  clear model_collection best_model;
  [test_frames_N, vars_N] = size(inferred_model.x_test);

  smart_chance = sum(inferred_model.x_test,1)./test_frames_N*100;
  %pred is vars by frames:
  PRED = CLL_1 > CLL_0;
  %GND is just inferred_model.x_test
  %x_test is frames by vars; how many correctly predicted as 1:
  TP = sum(PRED'.*inferred_model.x_test, 1)./sum(inferred_model.x_test,1)*100;
  TN = sum(~PRED'.*~inferred_model.x_test, 1)./sum(~inferred_model.x_test,1).*100;
  %accuracy = (tp+tn)/length(gnd);
  ACC(count,:) = (TP+TN)./vars_N;
  count = count + 1;
end
fprintf('saved: %s\n', fsave);
return

