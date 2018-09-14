% loopy data processing pipeline

rng(1000);
addpath(genpath('C:\Shuting\graph_ensemble\scripts'));
addpath('C:\Shuting\graph_ensemble\fwMatch-darpa\src\loopy_model');
addpath('C:\Shuting\Code\BCT\2016_01_16_BCT');
addpath('C:\Shuting\Code\CNMF DataProcessing\data_processing');

%% set parameters
param = struct();

param.expt_name = {'m21_d2_vis','m37_d2','pa_511510855_TF1','pa_511510650_TF1',...
     'pa_511507650_TF1','pa_511509529_TF1'};
param.ee = {{'all_high_add_neuron'},{'all_high_add_neuron'},...
    {'vis_high_add_neuron'},{'vis_high_add_neuron'},{'vis_high_add_neuron'},...
    {'vis_high_add_neuron'}};

% param.num_stim = [2,1];
param.ge_type = 'full'; % 'full', 'on', 'thresh'
param.savestr = ''; % 'add neuron'

% param.epth = 0.3; % quantile threshold of edge potential
param.epth = 0.05;

param.num_shuff = 100;
param.num_rand = 100;
param.k_seq = 2:6;
param.k = 3;
param.p = 0.05;

param.mc_minsz = 3;

param.data_path = 'C:\Shuting\graph_ensemble\data\';
param.shuff_path_base = [param.data_path 'shuffled\'];
param.result_path_base = 'C:\Shuting\graph_ensemble\results';
param.result_path.stats_path = [param.result_path_base '\stats\'];
param.fig_path.graph_prop = [param.result_path_base '\fig\graph_prop\'];
param.fig_path.core = [param.result_path_base '\fig\core\'];
param.fig_path.ens = [param.result_path_base '\fig\ens\'];
param.fig_path.pred = [param.result_path_base '\fig\pred\'];
param.fig_path.stats = [param.result_path_base '\fig\stats\'];
param.fig_path.opto_spont = [param.result_path_base '\fig\opto_spont\'];
param.fig_path.opto_stim = [param.result_path_base '\fig\opto_stim\'];
param.fig_path.pa = [param.result_path_base '\fig\pa\'];
param.fig_path.cc = [param.result_path_base '\fig\cc\'];
param.ccode_path = [param.result_path_base '\mycc.mat']; % color code
param.rwbmap = [param.result_path_base '\rwbmap.mat']; % red white blue map
param.graymap = [param.result_path_base '\graymap.mat']; % gray map
param.bluemap = [param.result_path_base '\bluemap.mat']; % blue map
param.redmap = [param.result_path_base '\redmap.mat']; % red map
param.redmap_light = [param.result_path_base '\redmap_light.mat']; % lighter red map
param.four_stim_cmap = [param.result_path_base '\four_stim_cmap.mat']; % four colors
param.bkrmap = [param.result_path_base '\bkrmap.mat']; % blue-black-red map

% threshold of high OSI cells
param.OSI_thresh = 0.7;

% parameters for reduced models
param.maxf = 1000; % max frames for reduced models
param.winsz = 100;

param.ndeg_bin_range = 0:0.01:0.25;
param.lcc_bin_range = 0:0.02:1;
param.cent_bin_range = 0:0.02:1;
param.mc_sz_bin_range = 0:1:15;
param.epsum_bin_range = -1:0.1:0.2;

% plotting parameters
param.linew = 0.5;

%% process data
% make shuffled data for CC graph
makeShuffledData(param);

% make cc graph
makeCCgraph(param);

% threshold edge potentials
threshCRFgraphs(param);

% threshold reduced models
for n = 1:length(param.expt_name)
    load([param.data_path param.expt_name{n} '\Pks_Frames.mat']);
    num_frame = length(Pks_Frame);
    trunc_vec = param.winsz:param.winsz:min([param.maxf,...
        floor(num_frame/param.winsz)*param.winsz]);
    for ii = 1:length(trunc_vec)
        param_rd = param;
        param_rd.expt_name = param.expt_name(n);
        param_rd.ee = {{[param.ee{n}{1} '_' num2str(trunc_vec(ii))]}};
        threshCRFgraphs(param_rd);
    end
end

% calculate cc and crf maximal cliques
calcGraphMC(param);

% calculate random graph maximal cliques
randGraphMC(param);

% % calculate cc and crf graph communities
% calcGraphComm(param);
% 
% % random graph properties for cc and crf
% randGraphStatsCC(param);
% randGraphStatsCRF(param);
% 
% % random graph comm and clique stats
% randGraphCommStats(param)

%% figures
% prediction with whole model
an_param = param;
an_param.expt_name = {'m21_d2_vis','m37_d2','pa_511510855_TF1','pa_511510650_TF1',...
     'pa_511507650_TF1','pa_511509529_TF1'};
an_param.ee = {{'all_high_add_neuron'},{'all_high_add_neuron'},...
    {'vis_high_add_neuron'},{'vis_high_add_neuron'},{'vis_high_add_neuron'},...
    {'vis_high_add_neuron'}};
fig2_crf_LL_pred_add_neuron(an_param);
% fig2_crf_LL_pred_add_neuron_multistim(an_param); % - not working well

% find ensembles
fig2_plot_ensemble_identification(param);

% compare CRF, SVD, OSI
figS2_SVD_OSI_pred(param);

% optimal ensembles
fig6_ensemble_reduction(param);

%% testing: finding ensembles in PA data
pa_param = param;
pa_param.expt_name = {'pa_511510855_TF1','pa_511510650_TF1',...
    'pa_511507650_TF1','pa_511507650_TF1'};
pa_param.ee = {{'vis_high_add_neuron'},{'vis_high_add_neuron'},...
    {'vis_high_add_neuron'}};

% fig2_plot_ensemble_identification(pa_param);

%% opto spont data
opto_param = param;
opto_param.ndeg_bin_range = 0.02:0.02:0.12;
opto_param.expt_name = {'m23_d1_opto_corrected'};
opto_param.ee = {{'high_pre_add_neuron','high_post_add_neuron'}};

fig8_opto_spont(opto_param);

%% opto stim data
opto_stim_param = param;
opto_stim_param.npot_bin_range = -1:0.05:1;
opto_stim_param.epot_bin_range = -1:0.05:1;
opto_stim_param.expt_name = {'m52_d1_opto'};
opto_stim_param.ee = {{'high_all_add_neuron'}};
opto_stim_param.ndeg_quant = 0.2;

fig7_opto_stim(opto_stim_param);

%% supplementary figures
% PA dataset
pa_param = param;
pa_param.expt_name = {'pa_511510855_1_TF1','pa_511510650_1_TF1',...
    'pa_511507650_1_TF1','pa_511509529_1_TF1','pa_511510670_1_TF1'};
pa_param.ee = {{'vis_high_add_neuron'},{'vis_high_add_neuron'},...
    {'vis_high_add_neuron'},{'vis_high_add_neuron'},{'vis_high_add_neuron'}};
pa_param.tf_seq = [1,2,4,8,15];
for n = 1:length(pa_param.expt_name)
    for m = 1:length(pa_param.tf_seq)
        pa_param.test{n}{m} = [pa_param.expt_name{n}(1:end-5) '23_TF' num2str(pa_param.tf_seq(m))...
             '_vis_high'];
    end
end
pa_param.qnoise = 0.3;
figS1_pa_whole_model_TF(pa_param);

crf_find_ensembles(pa_param);

% compare base model and hidden node model
cmp_param = param;
cmp_param.expt_name = {'m21_d2_vis','m37_d2','pa_511507650_1_TF1','pa_511509529_1_TF1',...
    'pa_511510650_1_TF1','pa_511510855_1_TF1'};
cmp_param.ee = {{'all_high','all_high_add_neuron'},{'all_high','all_high_add_neuron'},...
    {'vis_high','vis_high_add_neuron'},{'vis_high','vis_high_add_neuron'},...
    {'vis_high','vis_high_add_neuron'},{'vis_high','vis_high_add_neuron'}};
figS2_compare_add_neuron(cmp_param);

% compare number of neurons
% nn_param = param;
% nn_param.expt_name = {'m21_d2_vis'};
% nn_param.ee = {{'all_high','all_high_first_half_rand','all_high_second_half_rand'}};
% % nn_param.ee = {{'all_high','all_high_first_half_seq','all_high_second_half_seq'}};
% fig_cmp_neuron_num(nn_param);
nn_param = param;
nn_param.expt_name = {'m21_d2_vis'};
nn_param.ee = {{'all_high'}};
nn_param.rand_perc = 0.2:0.2:0.8;
gnReduceData(nn_param);
fig_cmp_neuron_num(nn_param);
