% loopy data processing pipeline

rng(1000);

%% set parameters
param = struct();
param.expt_name = {'m21_d2_vis','m37_d2'};
param.ee = {{'all_high_add_neuron'},{'all_high_add_neuron'}};
param.num_stim = [2,1];
param.cc_type = 'full';
param.comm_type = 'full';
param.ge_meth = 'full'; % 'full' or 'on'
param.savestr = 'add_neuron';

param.num_shuff = 100;
param.num_rand = 100;
param.k_seq = 2:6;
param.k = 3;
param.p = 0.05;

param.data_path = 'C:\Shuting\fwMatch\data\';
param.shuff_path_base = 'C:\Shuting\fwMatch\data\shuffled\';
param.result_path_base = 'C:\Shuting\fwMatch\results\';
param.result_path.stats_path = 'C:\Shuting\fwMatch\results\stats\';
param.fig_path.graph_prop = 'C:\Shuting\fwMatch\results\fig\graph_prop\';
param.fig_path.core = 'C:\Shuting\fwMatch\results\fig\core\';
param.fig_path.pred = 'C:\Shuting\fwMatch\results\fig\pred\';
param.fig_path.stats = 'C:\Shuting\fwMatch\results\fig\stats\';
param.ccode_path = 'C:\Shuting\fwMatch\results\mycc.mat'; % color code
param.rwbmap = 'C:\Shuting\fwMatch\results\rwbmap.mat'; % red white blue map

param.OSI_thresh = 0.4;

param.comm_sz_bin_range = 0:0.02:0.6;
param.comm_deg_bin_range = 0:0.1:5; %0:0.02:1;
param.comm_ov_bin_range = 0:0.02:0.6;
param.comm_mem_bin_range = 0:0.02:0.7;
param.ndeg_bin_range = 0:0.01:0.25;
param.lcc_bin_range = 0:0.02:1;
param.cent_bin_range = 0:0.02:1;
param.mc_sz_bin_range = 0:1:15;

param.linew = 0.5;

%% process data
% make shuffled data
makeShuffledData(param);

% make cc graph
makeCCgraph(param);

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

%%


