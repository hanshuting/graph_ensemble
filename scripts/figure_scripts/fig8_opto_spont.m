function [] = fig8_opto_spont(param)

% parameters
expt_name = param.expt_name;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.opto_spont;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
npot_bin_range = param.npot_bin_range;
epot_bin_range = param.epot_bin_range;
epsum_bin_range = param.epsum_bin_range;
linew = param.linew;

ndeg_bin_range = param.ndeg_bin_range;
lcc_bin_range = param.lcc_bin_range;
cent_bin_range = param.cent_bin_range;
mc_sz_bin_range = param.mc_sz_bin_range;

load(ccode_path);
graymap = load(param.graymap); % variable: cmap
bmap = load(param.bluemap); % variable: cmap

%% initialize
ndeg = cell(length(expt_name),2);
ndeg_in = cell(length(expt_name),2);
ndeg_cum = cell(length(expt_name),2);
ndeg_in_cum = cell(length(expt_name),2);

dens = cell(length(expt_name),2);
dens_in = cell(length(expt_name),2);

lcc = cell(length(expt_name),2);
lcc_in = cell(length(expt_name),2);
lcc_cum = cell(length(expt_name),2);
lcc_in_cum = cell(length(expt_name),2);

cent = cell(length(expt_name),2);
cent_in = cell(length(expt_name),2);
cent_cum = cell(length(expt_name),2);
cent_in_cum = cell(length(expt_name),2);

ep_sum = cell(length(expt_name),1);
ep_sum_cum = cell(length(expt_name),1);

npot_out = cell(length(expt_name),2);
npot_in = cell(length(expt_name),2);
npot_out_cum = cell(length(expt_name),2);
npot_in_cum = cell(length(expt_name),2);

epot_out = cell(length(expt_name),2);
epot_in = cell(length(expt_name),2);
epot_out_cum = cell(length(expt_name),2);
epot_in_cum = cell(length(expt_name),2);

% mc_num = zeros(length(expt_name),2);
% mc_sz = cell(length(expt_name),2);
% mc_sz_cum = cell(length(expt_name),2);

%% collect data from all experiments
for n = 1:length(expt_name)
    
    expt_ee = param.ee{n};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Stim_cells.mat']);
    num_node = size(Spikes,1);
    num_stim_cell = length(Stim_cells);
    nostim_cells = setdiff(1:num_node,Stim_cells);
    num_nostim_cell = num_node-num_stim_cell;
    
    pre_model = load([model_path expt_name{n} '_' expt_ee{1} '_loopy_best_model_' ge_type '.mat']);
    post_model = load([model_path expt_name{n} '_' expt_ee{2} '_loopy_best_model_' ge_type '.mat']);
    post_model.graph = full(post_model.graph);
    
    % convert to on edges
    pre_model.ep_on = getOnEdgePot(pre_model.graph,pre_model.G)';
    post_model.ep_on = getOnEdgePot(post_model.graph,post_model.G)';
    
    % extend coordinates for add neuron model
    coords = Coord_active;
    coords(end+1,:) = [0 max(coords(:,2))];
    coords(end+1,:) = [0 0];
    
    % plot pre and post models
    figure; set(gcf,'color','w','position',[2154 340 941 597])
    cc_range = [min(pre_model.ep_on(:)) max(pre_model.ep_on(:))];
    ep_range = [-1 0.2];
    subplot(1,2,1)
    plotGraphModelHighlightEP(pre_model.graph,coords,pre_model.ep_on,...
        cc_range,ep_range,graymap.cmap,[]);
    subplot(1,2,2)
    plotGraphModelHighlightEP(post_model.graph,coords,post_model.ep_on,...
        cc_range,ep_range,bmap.cmap,[]);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_ON_unnormalized_' ...
        ge_type '_pre_post_models.pdf']);

    % edge potential sum
    ep_sum{n,1} = sum(pre_model.ep_on,2);
    ep_sum{n,2} = sum(post_model.ep_on,2);
    ep_sum_cum{n,1} = calc_cum_dist(ep_sum{n,1},epsum_bin_range);
    ep_sum_cum{n,2} = calc_cum_dist(ep_sum{n,2},epsum_bin_range);
    
    % density
    dens{n,1} = sum(sum(pre_model.graph))/num_node/(num_node-1);
    dens{n,2} = sum(sum(post_model.graph))/num_node/(num_node-1);
    
    % node degree
    ndeg{n,1} = sum(pre_model.graph,2)/2/num_node;
    ndeg{n,2} = sum(post_model.graph,2)/2/num_node;
    ndeg_cum{n,1} = calc_cum_dist(ndeg{n,1},ndeg_bin_range);
    ndeg_cum{n,2} = calc_cum_dist(ndeg{n,2},ndeg_bin_range);
    
    % lcc
    lcc{n,1} = local_cluster_coeff(pre_model.graph);
    lcc{n,2} = local_cluster_coeff(post_model.graph);
    lcc_cum{n,1} = calc_cum_dist(lcc{n,1},lcc_bin_range);
    lcc_cum{n,2} = calc_cum_dist(lcc{n,2},lcc_bin_range);
    
    % centrality
    cent{n,1} = eigenvec_centrality(pre_model.graph);
    cent{n,2} = eigenvec_centrality(post_model.graph);
    cent_cum{n,1} = calc_cum_dist(cent{n,1},cent_bin_range);
    cent_cum{n,2} = calc_cum_dist(cent{n,2},cent_bin_range);

end

save([save_path 'opto_spont_prop.mat'],'-v7.3');

%% graph properties - nostim network
stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

figure; set(gcf,'color','w','position',[2006 320 691 425])

% edge pot sum
subplot(2,3,2)
plot_opto_spont_cum_hist(ep_sum_cum,mycc,epsum_bin_range,linew);
gcapos = get(gca,'position');
xlabel('sum(edge pot)');ylabel('log(p)');
set(gca,'position',gcapos);
legend off; box off

% density
boxwd = 0.2;
subplot(2,3,1);hold on;
dens_pre = cell2mat(dens(:,1)')*100;
h = boxplot(dens_pre,'positions',stepsz,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
dens_post = cell2mat(dens(:,2)')*100;
h = boxplot(dens_post,'positions',2*stepsz,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 3*stepsz]);
ylim([min([dens_pre,dens_post])-0.02 max([dens_pre,dens_post])]+0.02)
gcapos = get(gca,'position');
ylabel('density (%)')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);
box off

% node degree
subplot(2,3,3)
plot_opto_spont_cum_hist(ndeg_cum,mycc,ndeg_bin_range,linew);
gcapos = get(gca,'position');
xlabel('node degree');ylabel('log(p)');
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(2,3,4)
plot_opto_spont_cum_hist(lcc_cum,mycc,lcc_bin_range,linew);
gcapos = get(gca,'position');
xlabel('clustering coeff');ylabel('log(p)');
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(2,3,5)
plot_opto_spont_cum_hist(cent_cum,mycc,cent_bin_range,linew);
gcapos = get(gca,'position');
xlabel('centrality');ylabel('log(p)');
set(gca,'position',gcapos);
legend off; box off

saveas(gcf,[fig_path 'opto_spont_graph_prop_on_unnormalized_full.pdf']);


end