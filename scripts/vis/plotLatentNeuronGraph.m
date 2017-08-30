function [h] = plotLatentNeuronGraph(adjmat,coords,h)
% plot adjmat using the given node coordinates, with color specifying the
% weight (red: positive weight; blue: negative weight)

nodesz = 50;
r = 0.5; % line transparency
num_node = size(adjmat,1);

% add neuron
num_add = num_node-size(coords,1);
coords(end+1,:) = [0 max(coords(:,2))];
coords(end+1,:) = [0 0];
coords(end+1,:) = [max(coords(:,1)) 0];
coords(end+1,:) = [max(coords(:,1)) max(coords(:,2))];
coords = coords(1:num_node,:);

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

if nargin < 3
    h = figure;
    set(gcf,'color','w','position',[2056 296 396 344],'PaperPositionMode','auto');
end

hold on;
for i = 1:num_node
    
    % plot edge
    E = find(adjmat(i,:));
    
    if i == num_node-num_add+1
        cc = [1 0 0 r]; %red
    elseif i == num_node-num_add+2
        cc = [0 0 1 r]; %blue
    elseif i == num_node-num_add+3
        cc = [1 0 1 r]; %m
    elseif i == num_node-num_add+4
        cc = [0 1 1 r]; %c
    else
        cc = 0.7*[1 1 1 r];
    end
    
    if ~isempty(E)
        linkNode = coords(adjmat(i,:)~=0,:);
        crNode = repmat(coords(i,:),length(E),1);
        for j = 1:length(E)
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
        end
    end
    
end

for i = 1:num_node
    
    % plot nodes
    E = find(adjmat(i,:));
    
    if isempty(E)
        cc = 0.7*[1 1 1];
        mk = 'o';
    elseif any(E==num_node-num_add+1)
        cc = 'r';
        mk = 's';
    elseif any(E==num_node-num_add+2)
        cc = 'b';
        mk = 's';
    elseif any(E==num_node-num_add+3)
        cc = 'm';
        mk = 's';
    elseif any(E==num_node-num_add+4)
        cc = 'c';
        mk = 's';
    else
        cc = 0.7*[1 1 1];
        mk = 'o';
    end
    scatter(coords(i,1),coords(i,2),nodesz,mk,'markerfacecolor',cc,...
        'markeredgecolor','k','linewidth',1);
    
end

axis off equal tight

end
