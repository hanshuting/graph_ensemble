function [] = plotGraphModelHighlightSubset(adjmat,coords,indx,weightmat,cc_range)
% plot adjmat using the given node coordinates, with color specifying the
% weight (red: positive weight; blue: negative weight), size specifying
% node degrees
% INPUT:
%     adjmat: binary adjacency matrix
%     coords: N-by-2 matrix of node coordinates
%     weightmat: N-by-N matrix with edge weights
%     cc_range: color range, pass empty matrix [] if want to use default
%     values
%     cmap: 64x3 matrix of color map
% 
% Shuting Han, 2016

nodesz_max = 300;
nodesz_min = 10;
if isempty(cc_range)
    cc_range = max(abs(weightmat(:)))*1.01;
end
node_deg = sum(adjmat,2)/2;
nodesz_vec = node_deg/length(node_deg)*(nodesz_max-nodesz_min)+nodesz_min;

N = size(adjmat,1);

% set colormap
cmap = gray(64);
cmap_highlight = [zeros(64,1), (0:4:255)', 255*ones(64,1)]/255;

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

% plot edge
hold on;
for i = 1:N
    E = find(adjmat(i,:));
    if ~isempty(E)
        linkNode = coords(adjmat(i,:)~=0,:);
%         crNode = repmat(coords(i,:),length(E),1);
        for j = 1:length(E)
            cindx = ceil((weightmat(i,E(j))+cc_range)/(2*cc_range)*64);
            if cindx > 64; cindx = 64; end
            if cindx <=0; cindx = 1; end
            if ismember(i,indx) && ismember(j,indx)
                cc = cmap_highlight(cindx,:);
            else
                cc = cmap(cindx,:);
            end
            plot([coords(i,1),linkNode(j,1)]',[coords(i,2),linkNode(j,2)]','color',cc);
%             plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
%                 'color',cc);
        end
    end
end

colorbar;colormap(cmap);
caxis([-cc_range cc_range]);

% plot node
for i = 1:N
    if ismember(i,indx); cc = [0 0.5 1];
    else; cc = 0.7*[1 1 1]; end
    scatter(coords(i,1),coords(i,2),nodesz_vec(i),'markeredgecolor','k',...
        'markerfacecolor',cc,'linewidth',0.5);
end

axis off equal tight
hold off

end
