% assume "Coord_active" have been loaded
expt_name = 'mf2_v1';
ee = {'vis_01','vis_02','vis_03','vis_04','vis_all'};

for n = 1:length(ee)
    
%% load data
filename = [expt_name '_' ee{n}];
filepath = ['C:\Shuting\fwMatch\data\' expt_name '\'];
savepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
load([filepath filename '.mat']);

% calculate correlation
cc = corr(data);

% threshold correlation
cc_thresh = multithresh(cc(:));

% generate graph
graph = cc>cc_thresh;
graph = graph-diag(diag(graph));

% save result
save([savepath expt_name '_cc_model.mat'],'graph');

% 
% % plot graph
% plotGraphWithXY(graph,Coord_active,ones(size(graph)));
% 
% % community and centrality


end