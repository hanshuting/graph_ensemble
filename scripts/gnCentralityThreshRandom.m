% assume "Coord_active" have been loaded
expt_name = 'mf2_v1';
ee = {'vis_01','vis_02','vis_03','vis_04','vis_all'};
num_shuff = 100;

%% centrality
for n = 1:length(ee)
    
% load data
filename = [expt_name '_' ee{n} '_loopy_best_model'];
filepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
savepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
load([filepath filename '.mat']);

% number of node and edge
num_node = size(graph,1);
num_edge = sum(graph(:))/2;

% centrality from random graphs
node_cent_random = zeros(num_node,num_shuff);
for i = 1:num_shuff
    rand_graph = mkRandomGraph(num_node,num_edge);
    node_cent_random(:,i) = eigenvec_centrality(rand_graph);
end

% calculate single threshold
cent_thresh_single = zeros(num_node,1);
for i = 1:num_node
    ccdist = fitdist(node_cent_random(i,:)','normal');
    cent_thresh_single(i) = icdf(ccdist,0.95);
end

% calculate threshold all
ccdist = fitdist(node_cent_random(:),'normal');
cent_thresh_all = icdf(ccdist,0.95);

% save result
save([savepath exptname '_loopy_best_model_thresh_random.mat'],...
    'cent_thresh_all','cent_thresh_single','-v7.3');

end

%% community core
k = 3;
for n = 1:length(ee)
    
% load data
filename = [expt_name '_' ee{n} '_loopy_best_model'];
filepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
savepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
load([filepath filename '.mat']);

% number of node and edge
num_node = size(graph,1);
num_edge = sum(graph(:))/2;

% centrality from random graphs
comm_core_random = zeros(num_node,num_shuff);
for i = 1:num_shuff
    
%     rand_graph = mkRandomGraph(num_node,num_edge);
    rand_graph = linkRandomization(graph);
    
    % find community core cells
    [comm,~,~] = k_clique(k,rand_graph);
    cell_count = zeros(size(rand_graph,1),1);
    num_comm = length(comm);
%     for j = 1:num_comm
%         cell_count(comm{j}) = cell_count(comm{j})+1;
%     end
%     comm_core_random(:,i) = cell_count;
    comm_core_random(:,i) = histc(cell2mat(comm),1:num_node);
end

% calculate single threshold
comm_thresh_single = zeros(num_node,1);
for i = 1:num_node
    ccdist = fitdist(comm_core_random(i,:)','normal');
    comm_thresh_single(i) = icdf(ccdist,0.95);
end

% calculate threshold all
ccdist = fitdist(comm_core_random(:),'normal');
comm_thresh_all = icdf(ccdist,0.95);

% save result
save([savepath expt_name '_loopy_best_model_comm_thresh_random.mat'],...
    'comm_thresh_all','comm_thresh_single','-v7.3');

end
