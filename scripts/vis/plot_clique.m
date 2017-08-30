function [] = plot_clique(h,cliques,coords,nodesz)

set(h,'color','w');
scatter(coords(:,1),-coords(:,2),nodesz,'k','linewidth',1);
hold on;
scatter(coords(cliques,1),-coords(cliques,2),nodesz,'k',...
    'linewidth',1,'MarkerFaceColor','r');

axis off equal

end