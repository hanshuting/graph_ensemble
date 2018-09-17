function [] = plotHighlightSecondOrder(adjmat,coords,num_stim)
% highlight nodes and connecting edges specified by indx with specified color

load('C:\Shuting\graph_ensemble\results\mycc.mat');
ccstr = {mycc.red,mycc.blue,mycc.orange,mycc.purple};
% ccstr = {'r','b','m','c'};
num_node = size(adjmat,1);

% add neuron
num_add = num_node-size(coords,1);
coords(end+1,:) = [0 0];
coords(end+1,:) = [0 max(coords(:,2))];
coords(end+1,:) = [max(coords(:,1)) 0];
coords(end+1,:) = [max(coords(:,1)) max(coords(:,2))];
coords = coords(1:num_node,:);

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

nodesz = 50;
graycc = 0.7*[1 1 1];
num_node = size(adjmat,1);

figure

for n = 1:num_stim
    
    subplot(1,num_stim,n);
    hold on
    
    % plot edge
    for ii = 1:num_node-num_stim
        E = find(adjmat(ii,:));
        if ~isempty(E)
            linkNode = coords(adjmat(ii,:)~=0,:);
            crNode = repmat(coords(ii,:),length(E),1);
            for kk = 1:length(E)
                plot([crNode(kk,1),linkNode(kk,1)]',[crNode(kk,2),linkNode(kk,2)]',...
                    'color',graycc);
            end
        end
    end

    % highlight
    nodeSet = [];
    firstNode = find(adjmat(num_node-num_stim+n,:)~=0);
    nodeSet = unique([nodeSet,firstNode]);
    for jj = 1:length(firstNode)
        secondNode = find(adjmat(firstNode(jj),:)~=0);
        nodeSet = unique([nodeSet,secondNode]);
        if ~isempty(secondNode)
            linkNode = coords(secondNode,:);
            crNode = repmat(coords(firstNode,:),length(secondNode),1);
            for kk = 1:length(secondNode)
                plot([crNode(kk,1),linkNode(kk,1)]',[crNode(kk,2),linkNode(kk,2)]',...
                    'color',ccstr{n});
            end
        end
    end

    % plot node on top
    scatter(coords(:,1),coords(:,2),nodesz,'markerfacecolor',graycc,...
        'markeredgecolor','k','linewidth',1);

    % highlight
    scatter(coords(nodeSet,1),coords(nodeSet,2),nodesz,'markerfacecolor',...
        ccstr{n},'markeredgecolor','k','linewidth',1);

    axis off equal tight

end

end