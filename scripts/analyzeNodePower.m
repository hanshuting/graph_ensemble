% "node power"

EXPT = 'm21_sample_d2_d3';
EE = 'all_vis_01';
MODEL_NAME = 'local_model_4';
MODEL_TYPE = 'loopy';
DPATH = 'C:\Shuting\fwMatch\results\models';
SPATH = 'C:\Shuting\fwMatch\results\models';

%% load data
% raw data with coordinates
fload = sprintf('C:\Shuting\fwMatch\data\\%s.mat',EXPT);
load(fload);

% model data
fload = sprintf('%s\\%s\%s_%s_%s.mat',DPATH,EXPT,EXPT,EE,MODEL_NAME);
load(fload);

% shuffled data
fload = sprintf('%s\\shuffled\\shuffled_%s_%s_%s.mat',DPATH,EXPT,EE,MODEL_TYPE);
load(fload);

%%
% expt node power
node_power = abs(node_pot+sum(edge_pot,2)/2);

% shuffled threshold
num_node = length(node_pot);
num_shuffle = size(node_pot_shuffled,2);
node_power_shuffled = zeros(num_node,num_shuffle);
for i = 1:num_node
    for j = 1:num_shuffle
        node_power_shuffled(i,j) = abs(node_pot_shuffled(i,j)+...
            sum(squeeze(edge_pot_shuffled(i,:,j)))/2);
    end
end

% threshold
node_power_thresh = zeros(num_node,2);
for i = 1:num_node
    ccdist = fitdist(exp(node_power_shuffled(i,:)'),'normal');
    node_power_thresh(i,1) = icdf(ccdist,0.05);
    node_power_thresh(i,2) = icdf(ccdist,0.95);
end

node_power_thresh = log(node_power_thresh);
node_indx = node_power>node_power_thresh(:,2)|node_power<node_power_thresh(:,1);
node_indx = find(node_indx);

%% plot
h = figure;
set(h,'color','w');
scatter(Coord_active(:,1),-Coord_active(:,2),100,'k','linewidth',2);
hold on;
scatter(Coord_active(node_indx,1),-Coord_active(node_indx,2),100,'k',...
  'linewidth',2,'MarkerFaceColor','r');
axis off equal

%% eigenvector centrality
num_node = length(node_pot);
node_cent = eigenvec_centrality(graph);

% shuffle threshold
num_shuffle = size(node_pot_shuffled,2);
node_cent_shuffled = zeros(num_node,num_shuffle);
for i = 1:num_shuffle
    node_cent_shuffled(:,i) = eigenvec_centrality(logical(edge_pot_shuffled(:,:,i)));
end

% threshold
node_cent_thresh = zeros(num_node,1);
for i = 1:num_node
    ccdist = fitdist(node_cent_shuffled(i,:)','normal');
    node_cent_thresh(i) = icdf(ccdist,0.95);
end

node_indx = find(node_cent>node_cent_thresh);

%% plot centrality graph
cmap = jet(101);

h = figure;
set(h,'color','w');
hold on;
for i = 1:num_node
    cc = cmap(ceil((node_cent(i)-min(node_cent(:)))/...
        (max(node_cent)-min(node_cent))*100)+1,:);
    scatter(Coord_active(i,1),-Coord_active(i,2),100,'k','linewidth',2,...
        'MarkerFaceColor',cc);
end
axis off equal
colorbar;colormap(jet);
