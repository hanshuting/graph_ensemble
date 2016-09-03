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
  MODEL_TYPE = '_loopy';
  %ALL_EXPT = {'0', '45', '90', '135', '180', '225'};
  DEPTH = 1;
  ORIENT = 90;
  EXPT = sprintf('steph_%d_%02d%s', ORIENT, DEPTH, MODEL_TYPE); %'steph_135_01_tree';
  if(~isempty(who('SCRIPT')))
    warning(['You are calling from a script but not all vars are ' ...
             'set!!!']);
  end
end


%% extract tuning orientation from the python csv files:
EE = {'482923718', '482924833', '483020038', ...
      '483020476', '483056972', '483059231', ...
      '483061828'};

if(~exist('~/data/CAM/steph_all_tunings.mat', 'file'))
  c = 1;
  for e = EE
    expt = e{1};
    fload = sprintf('~/data/steph_225/tuning_%s_mean_preStim_responses.csv',expt);
    fprintf('reading: %s\n', fload);
    %rois_N by 3: index, orient, frequency
    data = dlmread(fload);
    TUNES{c} = data(:,2:3);
    c = c + 1;
    clear data;
  end
  
  fsave = sprintf('~/data/CAM/steph_all_tunings.mat');
  save(fsave, 'TUNES');
else
  load('~/data/CAM/steph_all_tunings.mat');
end

%% ROIS:
load('~/data/CAM/steph_all_rois.mat', 'ALL_ROIS', 'EE', ...
     'DEPTH_LABELS');


%% reliable neurons:
load('~/data/steph_225_tree/reliable_neurons.mat', 'RELIABLE');


%% models:

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


positions = zeros(vars_N,2);
for v = 1:vars_N
  variable_names{v} = sprintf('N%d',v);
  positions(v,:) = ALL_ROIS{DEPTH}(v,:);
end

positions(:,2) = 512 - positions(:,2);
SPATH = sprintf('~/data/steph_225_tree/%s_neuron_analysis', EXPT);
if(~exist(SPATH, 'dir'))
  mkdir(SPATH);
end
fsave = [SPATH '/bm_adj_' EXPT '.csv'];
csvwrite(fsave, full(inferred_model.structure));
fprintf('saved: %s\n', fsave);

%edges only for reliable neurons:
rel_adj = full(inferred_model.structure);
zero_out_nodes = setdiff(1:vars_N, RELIABLE{DEPTH});
rel_adj(:,zero_out_nodes) = 0;
fsave = [SPATH '/rel_adj_' EXPT '.csv'];
csvwrite(fsave, rel_adj);
fprintf('saved: %s\n', fsave);


fsave = [SPATH '/rois_xy_' EXPT '.csv'];
csvwrite(fsave, positions);
fprintf('saved: %s\n', fsave);

%reliable = b, reliable + tuned = g, rel + sampled = o, rel + tun +
%sampled = r
colors = ones(vars_N,3);
TUNED_EXPT = find(TUNES{DEPTH}(:,1) == ORIENT);
%SAMPLED_EXPT = find(ave_X >= quantile(ave_X, 0.9));
%doing one depth at a time
for i = 1:vars_N
  rel = ismember(i, RELIABLE{DEPTH});
  tuned = ismember(i, TUNED_EXPT);
  sampled = false; %ismember(i, SAMPLED_EXPT);
  if( rel & tuned & sampled)
    colors(i,:) = [1 0 0];
  elseif(rel & sampled )
    colors(i,:) = [1 0.6471 0];
  elseif(rel & tuned)
    colors(i,:) = [0 1 0];
  elseif(rel)
    colors(i,:) = [0 0 1];
  end
end

fsave = [SPATH '/colors_' EXPT '.csv'];
csvwrite(fsave, colors);
fprintf('saved: %s\n', fsave);



