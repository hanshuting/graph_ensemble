
% assume "graph", "edge_pot" and "Coord_active" have been loaded
expt_name = 'm21_d2_vis';
% ee = {'vis_01_all','vis_02_all','vis_03_all','vis_04_all'};
ee = {'01'};

for ii = 1:length(ee)
    
%% load data
filename = [expt_name '_' ee{ii} '_loopy_best_model'];
filepath = ['C:\Shuting\fwMatch\results\' expt_name '\models\'];
savepath = ['C:\Shuting\fwMatch\results\' expt_name '\fig\'];

load([filepath filename '.mat']);
M = full(graph);
num_node = size(graph,1);

%% plot best model
h = plotGraphWithXY(graph,Coord_active,edge_pot);
saveas(h,[savepath filename '_graph']);

%% find cliques and communities
k = 3;
[comm,cliques,~] = k_clique(k,M);

% find the cells that appear in multiple communities
cell_count = histc(cell2mat(comm),1:num_node);
core_indx = find(cell_count>1);

% plot
h = figure;
set(h,'color','w');
nodesz = 200;
scatter(Coord_active(:,1),-Coord_active(:,2),nodesz,'k','linewidth',2);
hold on;
scatter(Coord_active(core_indx,1),-Coord_active(core_indx,2),nodesz,'k',...
  'linewidth',2,'MarkerFaceColor','r');
for i = 1:length(core_indx)
    text(Coord_active(core_indx(i),1)-5,-Coord_active(core_indx(i),2),...
        num2str(cell_count(core_indx(i))),'color','y','fontsize',10,'fontweight','bold');
end
axis off equal
title('core cells');

saveas(h,[savepath filename '_k_' num2str(k) '_community_core']);

%% thresholded community membership
thresh = 8;

comm_sz = cellfun('length',comm);
cell_count = histc(cell2mat(comm(comm_sz>=thresh)),1:num_node);
core_indx = find(cell_count>1);

% plot
h = figure;
set(h,'color','w');
nodesz = 200;
scatter(Coord_active(:,1),-Coord_active(:,2),nodesz,'k','linewidth',2);
hold on;
scatter(Coord_active(core_indx,1),-Coord_active(core_indx,2),nodesz,'k',...
  'linewidth',2,'MarkerFaceColor','r');
for i = 1:length(core_indx)
    text(Coord_active(core_indx(i),1)-5,-Coord_active(core_indx(i),2),...
        num2str(cell_count(core_indx(i))),'color','y','fontsize',10,'fontweight','bold');
end
axis off equal
title('core cells');

saveas(h,[savepath filename '_k_' num2str(k) '_community_core_thresh']);

%% centrality analysis
node_cent = eigenvec_centrality(graph);

% plot centrality graph
num_node = length(node_pot);
cmap = jet(101);
h = figure;
set(h,'color','w');
hold on;
for i = 1:num_node
    cc = cmap(ceil((node_cent(i)-min(node_cent(:)))/...
        (max(node_cent)-min(node_cent))*100)+1,:);
    scatter(Coord_active(i,1),-Coord_active(i,2),nodesz,'k','linewidth',2,...
        'MarkerFaceColor',cc);
end
axis off equal
colorbar;colormap(jet);
title('centrality')

saveas(h,[savepath filename '_centrality']);

end