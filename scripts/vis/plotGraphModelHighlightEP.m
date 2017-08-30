function [] = plotGraphModelHighlightEP(adjmat,coords,weightmat,cc_range,...
    ep_range,edge_map,node_map)
% plot adjmat using the given node coordinates, with color specifying the
% weight (red: positive weight; blue: negative weight), size specifying
% node degrees
% INPUT:
%     adjmat: binary adjacency matrix
%     coords: N-by-2 matrix of node coordinates
%     weightmat: N-by-N matrix with edge weights
%     cc_range: color range, pass empty matrix [] if want to use default
%     values
% 
% Shuting Han, 2016

if isempty(edge_map)
%     cmap = load('C:\Shuting\graph_ensemble\results\rbmap.mat');
%     cmap = struct2array(cmap);
    edge_map = jet(64);
end
if isempty(node_map)
    node_map = jet(64);
end

nodesz_max = 300;
nodesz_min = 10;
if isempty(cc_range)
    cc_range = [min((weightmat(:))) max((weightmat(:)))];
end
node_deg = sum(adjmat,2)/2;
% nodesz_vec = (node_deg-min(node_deg))./(max(node_deg)-min(node_deg)).*...
%     (nodesz_max-nodesz_min)+nodesz_min;
nodesz_vec = node_deg/length(node_deg)*(nodesz_max-nodesz_min)+nodesz_min;

gcapos = get(gca,'position');
axis off

N = size(adjmat,1);

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

% plot edge
ax1 = axes;
for i = 1:N
    E = find(adjmat(i,:));
    if ~isempty(E)
        linkNode = coords(adjmat(i,:)~=0,:);
        crNode = repmat(coords(i,:),length(E),1);
        for j = 1:length(E)
            cindx = ceil((weightmat(i,E(j))-cc_range(1))/(cc_range(2)-cc_range(1))*64);
            if cindx > 64
                cindx = 64;
            end
            if cindx <=0
                cindx = 1;
            end
            cc = edge_map(cindx,:);
            plot(ax1,[crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
            hold on;
        end
    end
end

axis off equal tight
colormap(ax1,edge_map);colorbar(ax1,'eastoutside');
caxis([cc_range(1) cc_range(2)]);
set(ax1,'position',gcapos);

% sum up edge potentials
eps = nansum(weightmat,2);%./node_deg;
if isempty(ep_range)
    ep_range = [min((eps(:))) max((eps(:)))];
end

% plot node
ax2 = axes;
hold on;
for i = 1:N
    cindx = ceil((eps(i)-ep_range(1))/(ep_range(2)-ep_range(1))*64);
    if cindx<=0 || isnan(cindx)
        cindx = 1;
    elseif cindx >= 64
        cindx = 64;
    end
    if sum(adjmat(i,:))~=0
        scatter(ax2,coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
            'markerfacecolor',node_map(cindx,:),'linewidth',0.5);
    else
        scatter(ax2,coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
            'markerfacecolor','w','linewidth',0.5);
    end
end

axis off equal tight
hold off

colormap(ax2,node_map);colorbar(ax2,'southoutside');
caxis([ep_range(1) ep_range(2)]);
set(ax2,'position',gcapos);

linkaxes([ax1,ax2]);

end
