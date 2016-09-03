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
DEPTH = 2;


models_N = length(ALL_EXPT);
LLs = zeros(models_N);

test_samples = [];
gnd = [];

for e = 1:length(ALL_EXPT)
  EXPT = ALL_EXPT{e};
  %we don't need the best model; all the models use the same
  %train/test split.
  %fload = sprintf('%s/%s%s%s/results/model_collection.mat', LPATH, ...
  fload = sprintf('%s/%s%s_%02d%s/results/result1.mat', LPATH, ...
                  MODEL_PREFIX, EXPT, DEPTH, MODEL_TYPE);
  load(fload);
  fprintf('loaded model %s\n', fload);

  test_samples = [test_samples; model_collection.x_test];
  gnd = [gnd, str2num(ALL_EXPT{e}).*ones(1,size(model_collection.x_test,1))];
  clear model_collection;

end

fsave = sprintf('%s/test_samples_d%02d.mat', SPATH, DEPTH);
save(fsave, 'test_samples', 'gnd');
fprintf('SAVED: %s\n', fsave);