% dbstop error;

% Diary
%[~, name] = system('hostname');
%% strip the newline
%name(end) = [];
%diaryFile = [name ' ' datestr(now()) '.diary']
%diary(diaryFile);

basepath = pwd;

addpath(fullfile(basepath,'src'))
addpath(fullfile(basepath,'thirdparty'))
addpath(fullfile(basepath,'thirdparty','QPBO-v1.32.src'))
addpath(fullfile(basepath,'src','overcomplete'))

% Experiment Framework
addpath(fullfile(basepath,'src','expt_framework'))

% Structure Learning code
addpath(fullfile(basepath,'src','structure_learning'))
addpath(fullfile(basepath,'src','loopy_model'))
addpath(fullfile(basepath,'src','loopy_model','structure_learning'))
addpath(fullfile(basepath,'src','loopy_model','parameter_estimation'))

% Paramter Estimation
addpath(fullfile(basepath,'src','paramter_estimation'))

% ensemble tools
addpath(fullfile(basepath,'src','util'))
addpath(fullfile(basepath,'src','graphs'))
addpath(fullfile(basepath,'src','core'))
addpath(fullfile(basepath,'src','vis'))

display(sprintf('PATH is now set!'));
