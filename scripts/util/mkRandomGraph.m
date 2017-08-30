function [graph] = mkRandomGraph(N,E)
% make Erdos-Renyi random graphs
% N: number of nodes; E: number of edges
% 
% Shuting Han, 2016

% generate full edge list
[xx,yy] = meshgrid(1:N,1:N);
xx = xx(:);
yy = yy(:);
uindx = xx<yy;
edge_list = [xx(uindx),yy(uindx)];

% randomly pick E edges
edge_indx = randperm(size(edge_list,1),E);
edge_list = edge_list(edge_indx(:),:);

% put back to graph
graph = false(N,N);
eindx = sub2ind([N,N],edge_list(:,1),edge_list(:,2));
graph(eindx) = true;

% make it a full matrix
graph = graph+graph';

end