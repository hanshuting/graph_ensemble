function [h] = plotLatentNeuronMC(adjmat,coords)
% plot adjmat using the given node coordinates, with color specifying the
% weight (red: positive weight; blue: negative weight)

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

% find maximal cliques
mc = maximalCliques(adjmat);
mc_nz = cell(size(mc,2),1);
for n = 1:size(mc,2)
    mc_nz{n} = find(mc(:,n));
end

% find mc connected with latent nodes
mc_ln = cell(num_add,1);
indx = false(num_add,size(mc,2));
for n = 1:num_add
    indx(n,:) = logical(mc(num_node-num_add+n,:)~=0);
    mc_ln{n} = mc_nz(indx(n,:));
end

% plot
h = figure;
set(gcf,'color','w','position',[2056 459 248 587],'PaperPositionMode','auto');
for n = 1:num_add
    subplot(num_add,1,n)
    plotHighlightWithEdge(adjmat,coords,mc_ln{n},ccstr{n});
end

end
