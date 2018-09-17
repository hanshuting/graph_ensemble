function [] = plotHighlightWithEdge(adjmat,coords,mc_ln,cc)
% highlight nodes and connecting edges specified by indx with specified color

nodesz = 50;
r = 0.5; % line transparency
graycc = 0.7*[1 1 1];
num_node = size(adjmat,1);

hold on

%% plot edge
for ii = 1:num_node
    E = find(adjmat(ii,:));
    if ~isempty(E)
        linkNode = coords(adjmat(ii,:)~=0,:);
        crNode = repmat(coords(ii,:),length(E),1);
        for j = 1:length(E)
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',graycc);
        end
    end
end

% highlight
for ii = 1:length(mc_ln)
    cr_mc = mc_ln{ii};
    for jj = 1:length(cr_mc)
        E = setdiff(cr_mc,jj);
        linkNode = coords(E,:);
        crNode = repmat(coords(cr_mc(jj),:),length(E),1);
        for j = 1:length(E)
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
        end
    end
end

%% plot node on top
for ii = 1:num_node
    scatter(coords(ii,1),coords(ii,2),nodesz,'markerfacecolor',graycc,...
        'markeredgecolor','k','linewidth',1);
end

% highlight
for ii = 1:length(mc_ln)
    scatter(coords(mc_ln{ii},1),coords(mc_ln{ii},2),nodesz,'markerfacecolor',...
        cc,'markeredgecolor','k','linewidth',1);
end

axis off equal tight

end