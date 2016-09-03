
EXPT = 'm21_sample_d2_d3';
EE = 'all_vis_02';
MODEL_NAME = 'local_model_2';
MODEL_TYPE = 'loopy';
DPATH = 'C:\Shuting\fwMatch\results\models';
SPATH = 'C:\Shuting\fwMatch\results\models';

% load model
fload = sprintf('%s\\%s\\%s_%s_%s.mat',DPATH,EXPT,EXPT,EE,MODEL_NAME);
load(fload);
fload = sprintf('%s\\shuffled\\shuffled_%s_%s_%s.mat',DPATH,EXPT,EE,MODEL_TYPE);
load(fload);

% threshold
edge_pot = edge_pot.*(edge_pot<edge_thresh(:,:,1)|edge_pot>edge_thresh(:,:,2));
graph = edge_pot~=0;

% save
fsave = sprintf('%s\\thresh\\thresh_%s_%s_%s.mat',SPATH,EXPT,EE,MODEL_NAME);
save(fsave,'edge_pot','graph','-v7.3');

