function [] = visGraphCirc(adjmat,varargin)
% Visualize graph adjacency matrix with a circle layout in current figure
% Options can be specified with name value pairs:
%     'nodeColors': ['c'] A cell array or n-by-3 matrix specifying colors
%                   for the nodes. 
%     'edgeColors': An n-by-3 cell array listing {fromNode,toNode,color}
%                   for each row. You can list only the n < numel(edges) 
%                   edges you want to color. You can specify the text 'all'
%                   in place of toNode to mean all nodes, i.e. {fromNode,
%                   'all',color}. You can specify the color to be NaN in
%                   order to not plot a specific edge.
% Call the function without name value pairs to use default color settings.
% 
% Shuting Han, 2016

%% parse input
num_node = size(adjmat,1);
edge_list = zeros(sum(adjmat(:))/2,2);
[edge_list(:,2),edge_list(:,1)] = find(tril(adjmat));

ip = inputParser;
ip.StructExpand = true;
ip.KeepUnmatched = true;
ip.addParameter('ifvisualize',true,@islogical);
ip.addParameter('nodeColors',repmat([1 1 0.6],size(adjmat,1),1));
ip.addParameter('edgeColors',[num2cell(edge_list),mat2cell(repmat([0 0 0],...
    sum(adjmat(:))/2,1),ones(size(edge_list,1),1),3)]);
ip.parse(varargin{:});

if size(adjmat,1)~=size(adjmat,2)
    error('adjacency matrix must be symmetric')
end

params = ip.Results;

%% get node position on circle
boxsz = 110;
rr = 100;
node_cent = zeros(num_node,2);
for n = 1:num_node
    node_cent(n,:) = rr*[sin(n/num_node*2*pi),cos(n/num_node*2*pi)];
end

%% plot
maxnodesz = 500;
minnodesz = 10;
linew = 0.5;

hold on

% plot edges
for n = 1:size(edge_list,1)
    node1 = edge_list(n,1);
    node2 = edge_list(n,2);
    edge_indx = cellfun(@(x) x==node1,params.edgeColors(:,1)) & ...
        (cellfun(@(x) sum(x==node2,2),params.edgeColors(:,2))|...
        cellfun(@(x) strcmp(x,'all'),params.edgeColors(:,2)));
    edge_indx = edge_indx | cellfun(@(x) x==node2,params.edgeColors(:,1)) & ...
        (cellfun(@(x) sum(x==node1,2),params.edgeColors(:,2))|...
        cellfun(@(x) strcmp(x,'all'),params.edgeColors(:,2)));
    if any(edge_indx~=0)
        cc = params.edgeColors{edge_indx,3};
        cc = cc(1,:);
        if ~isnan(cc)
            plot([node_cent(node1,1),node_cent(node2,1)],[node_cent(node1,2),...
                node_cent(node2,2)],'color',cc,'linewidth',linew);
        end
    end
end

% plot nodes
for n = 1:num_node
    scatter(node_cent(n,1),node_cent(n,2),sum(adjmat(n,:))/num_node*maxnodesz+minnodesz,...
        params.nodeColors(n,:),'filled','MarkerEdgeColor','k','LineWidth',linew);
end

xlim([-boxsz boxsz]); ylim([-boxsz boxsz]);
axis equal tight off

end