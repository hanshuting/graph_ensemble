dbstop if error;

% Diary
%[~, name] = system('hostname');
%% strip the newline
%name(end) = [];
%diaryFile = [name ' ' datestr(now()) '.diary']
%diary(diaryFile);

addpath src/;
addpath run/;
addpath thirdparty/;
addpath thirdparty/blossom;
addpath thirdparty/pmtk3-master; startup_pmtk3;
%addpath thirdparty/csa++;
addpath thirdparty/permanent;
addpath thirdparty/bp_permanent;
addpath thirdparty/cm_and_cb_utilities;
addpath(genpath('thirdparty/JustinsGraphicalModelsToolboxPublic'));
addpath thirdparty/QPBO-v1.32.src;
addpath thirdparty/metis-5.1.0/metismex;
addpath thirdparty/fastAUC/;
%addpath src/fwBipartite;
addpath src/minReg;
addpath src/overcomplete;

% Commented because it is causing warnings for everyone
%addpath /proj/learning/software/lightspeed;
%addpath /proj/learning/software/matlab_xunit/xunit;
%addpath /Users/kuitang/Documents/MATLAB/toolboxes/lightspeed;

% Experiment Framework (added by Henrique)
addpath src/expt_framework

% Structure Learning code (added by Henrique)
addpath src/structure_learning
addpath src/structure_learning/simple_infer_structure/  
addpath src/gibbs_sampler/
addpath src/simple_graph_visualization/
addpath src/loopy_model/
addpath src/loopy_model/structure_learning
addpath src/loopy_model/parameter_estimation
addpath src/baseline_models/

% Paramter Estimation (added by Rashmi)
addpath src/paramter_estimation/;
addpath src/paramter_estimation/Drawing/;
addpath thirdparty/libdai/matlab/;

% Data Parsing (added by Liang)
addpath src/DataParsing/;

% Financial data (added by Liao)
addpath src/FinMatlabFunctions;

% javaaddpath thirdparty/lingpipe-4.1.0.jar;

display(sprintf('PATH is now set!'));
