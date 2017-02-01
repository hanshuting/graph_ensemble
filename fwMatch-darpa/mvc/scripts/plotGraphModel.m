function [] = plotGraphModel(adjmat,coords,weightmat,cc_range)
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

rbmap = load('C:\Shuting\fwMatch\results\rbmap.mat');
rbmap = struct2array(rbmap);
nodesz_max = 100;
nodesz_min = 10;
if isempty(cc_range)
    cc_range = max(abs(weightmat(:)))*1.01;
end
node_deg = sum(adjmat,2)/2;
nodesz_vec = (node_deg-min(node_deg))./(max(node_deg)-min(node_deg)).*...
    (nodesz_max-nodesz_min)+nodesz_min;

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
            cc = rbmap(ceil((weightmat(i,E(j))+cc_range)/(2*cc_range)*64),:);
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
            hold on;
        end
    end
end

colorbar;colormap(rbmap);
caxis([-cc_range cc_range]);

% plot node
for i = 1:N
    scatter(coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
        'markerfacecolor',[0.7 0.7 0.7],'linewidth',1.5);
end

axis off
hold off

end
