function [] = plotCoreOverlay(coords,indx1,indx2,cc1,cc2,rr)
% plot scatter plot neurons, highlight with coloring

nodesz = 50;
if nargin < 6
    rr = 1;
end

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

% find shared nodes
sind = intersect(indx1,indx2);

% plot
hold on;
% contour
scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',0.5);
% core 1
scatter(coords(indx1,1),coords(indx1,2),nodesz,'markeredgecolor','k',...
    'markerfacecolor',cc1,'linewidth',0.5,'markerfacealpha',rr);
% core 2 unique
scatter(coords(setdiff(indx2,sind),1),coords(setdiff(indx2,sind),2),nodesz,...
    'markeredgecolor','k','markerfacecolor',cc2,'linewidth',0.5,'markerfacealpha',rr);
% shared contour
scatter(coords(sind,1),coords(sind,2),nodesz,'markeredgecolor',cc2,...
    'markerfacecolor','w','linewidth',2.5,'markerfacealpha',0);

axis equal tight off

end
