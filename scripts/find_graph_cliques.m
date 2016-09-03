
% assume "graph", "edge_pot" and "Coord_active" have been loaded
expt_name = 'mf2_v1';
ee = 'vis_01';
filename = [expt_name '_' ee '_loopy_best_model'];
filepath = ['C:\Shuting\fwMatch\results\models\' expt_name '\'];
savepath = ['C:\Shuting\fwMatch\results\fig\' expt_name '\'];

load([filepath filename '.mat']);
M = full(graph);

%% find cliques and communities
k = 3;
[comm,cliques,~] = k_clique(k,M);

% plot cliques
for i = 1:length(cliques)
    h = figure;
    set(h,'color','w');
    scatter(Coord_active(:,1),-Coord_active(:,2),100,'k','linewidth',2);
    hold on;
    scatter(Coord_active(cliques{i},1),-Coord_active(cliques{i},2),100,'k',...
      'linewidth',2,'MarkerFaceColor','r');
    axis off equal
    title(['clique' num2str(i)]);
    saveas(h,[savepath filename '/' filename '_clique_' num2str(i)]);
end

% plot communities
for i = 1:length(comm)
    h = figure;
    set(h,'color','w');
    scatter(Coord_active(:,1),-Coord_active(:,2),100,'k','linewidth',2);
    hold on;
    scatter(Coord_active(comm{i},1),-Coord_active(comm{i},2),100,'k',...
      'linewidth',2,'MarkerFaceColor','r');
    axis off equal
    title(['community' num2str(i)]);
    saveas(h,[savepath filename '/' filename '_k_' num2str(k) '_community_' num2str(i)]);
end

% plot in a single graph
% cliques
num_cliques = length(cliques);
N = ceil(sqrt(num_cliques));
M = ceil(num_cliques/N);
figure;set(gcf,'color','w')
for i = 1:num_cliques
    h = subplot(M,N,i);
    scatter(h,Coord_active(:,1),-Coord_active(:,2),10,'k','linewidth',0.5);
    hold on;
    scatter(h,Coord_active(cliques{i},1),-Coord_active(cliques{i},2),10,'k',...
      'linewidth',0.5,'MarkerFaceColor','r');
    axis off equal
    title(num2str(i));
end
saveas(gcf,[savepath filename '/' filename '_clique_all']);

% communities
num_comm = length(comm);
N = ceil(sqrt(num_comm));
M = ceil(num_comm/N);
figure;set(gcf,'color','w')
for i = 1:num_comm
    h = subplot(M,N,i);
    scatter(h,Coord_active(:,1),-Coord_active(:,2),10,'k','linewidth',0.5);
    hold on;
    scatter(h,Coord_active(comm{i},1),-Coord_active(comm{i},2),10,'k',...
      'linewidth',0.5,'MarkerFaceColor','r');
    axis off equal
    title(num2str(i));
end
saveas(gcf,[savepath filename '/' filename '_k_' num2str(k) '_community_all']);

%% find the cells that appear in multiple communities
cell_count = zeros(size(graph,1),1);
num_comm = length(comm);
for i = 1:num_comm
    cell_count(comm{i}) = cell_count(comm{i})+1;
end
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
saveas(h,[savepath filename '_core_cells']);


