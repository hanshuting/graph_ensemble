function [] = plotGraphModelHighlight(adjmat,coords,weightmat,indx,cc_range,cmap)
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

if nargin < 6
    cmap = load('C:\Shuting\graph_ensemble\results\rbmap.mat');
    cmap = struct2array(cmap);
    cmap = jet(64);
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

N = size(adjmat,1);

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

% plot edge
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
            cc = cmap(cindx,:);
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
            hold on;
        end
    end
end

colorbar;colormap(cmap);
caxis([cc_range(1) cc_range(2)]);

% plot node
for i = 1:N
    if sum(indx==i)~=0
        scatter(coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
            'markerfacecolor','r','linewidth',0.5);
    else
        scatter(coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
            'markerfacecolor',[0.7 0.7 0.7],'linewidth',0.5);
    end
end

axis off equal tight
hold off

end
