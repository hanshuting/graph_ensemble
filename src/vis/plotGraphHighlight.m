function [] = plotGraphHighlight(coords,indx,colorstr,rr)
% plot scatter plot neurons, highlight with coloring

if nargin < 4
    rr = 1;
end

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

nodesz = 30;
scatter(coords(indx,1),coords(indx,2),nodesz,'markeredgecolor','k',...
    'markerfacecolor',colorstr,'linewidth',0.5,'markerfacealpha',rr);
hold on;
scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',0.5);
axis equal tight off
% hold off

end
