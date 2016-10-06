function [] = fig7_opto_spont(param)

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
linew = param.linew;

ndeg_bin_range = param.ndeg_bin_range;
lcc_bin_range = param.lcc_bin_range;
cent_bin_range = param.cent_bin_range;
mc_sz_bin_range = param.mc_sz_bin_range;

load(ccode_path);

%% initialize
ndeg_all = cell(length(expt_name),2);
ndeg_in = cell(length(expt_name),2);
ndeg_all_cum = cell(length(expt_name),2);
ndeg_in_cum = cell(length(expt_name),2);

dens_all = cell(length(expt_name),2);
dens_in = cell(length(expt_name),2);

lcc_all = cell(length(expt_name),2);
lcc_in = cell(length(expt_name),2);
lcc_all_cum = cell(length(expt_name),2);
lcc_in_cum = cell(length(expt_name),2);

cent_all = cell(length(expt_name),2);
cent_in = cell(length(expt_name),2);
cent_all_cum = cell(length(expt_name),2);
cent_in_cum = cell(length(expt_name),2);

npot_all = cell(length(expt_name),2);
npot_in = cell(length(expt_name),2);
npot_all_cum = cell(length(expt_name),2);
npot_in_cum = cell(length(expt_name),2);

epot_all = cell(length(expt_name),2);
epot_in = cell(length(expt_name),2);
epot_all_cum = cell(length(expt_name),2);
epot_in_cum = cell(length(expt_name),2);

mc_num = zeros(length(expt_name),2);
mc_sz = cell(length(expt_name),2);
mc_sz_cum = cell(length(expt_name),2);

%% collect data from all experiments
for n = 1:length(expt_name)
    
    expt_ee = param.ee{n};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Stim_cells.mat']);
    num_node = size(Spikes,1);
    num_stim_cell = length(Stim_cells);
    
    pre_model = load([model_path expt_name{n} '_' expt_ee{1} '_loopy_best_model_' ge_type '.mat']);
    post_model = load([model_path expt_name{n} '_' expt_ee{2} '_loopy_best_model_' ge_type '.mat']);
    
    % plot pre and post models
    figure; set(gcf,'color','w','position',[2154,119,1006,818])
    subplot(2,2,1);
    plotGraphModel(pre_model.graph,Coord_active,pre_model.edge_pot,[]);
    subplot(2,2,2);
    plotGraphModel(post_model.graph,Coord_active,post_model.edge_pot,[]);
    subplot(2,2,3);
    plotGraphModel(pre_model.graph(Stim_cells,Stim_cells),Coord_active...
        (Stim_cells,:),pre_model.edge_pot(Stim_cells,Stim_cells),[]);
    subplot(2,2,4);
    plotGraphModel(post_model.graph(Stim_cells,Stim_cells),Coord_active...
        (Stim_cells,:),post_model.edge_pot(Stim_cells,Stim_cells),[]);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
        ge_type '_pre_post_models.pdf']);
    
    % collect node and edge potentials
    npot_all{n,1} = pre_model.node_pot;
    npot_all{n,2} = post_model.node_pot;
    npot_in{n,1} = pre_model.node_pot(Stim_cells);
    npot_in{n,2} = post_model.node_pot(Stim_cells);
    npot_all_cum{n,1} = calc_cum_dist(npot_all{n,1},npot_bin_range);
    npot_all_cum{n,2} = calc_cum_dist(npot_all{n,2},npot_bin_range);
    npot_in_cum{n,1} = calc_cum_dist(npot_in{n,1},npot_bin_range);
    npot_in_cum{n,2} = calc_cum_dist(npot_in{n,2},npot_bin_range);
    
    epot_all{n,1} = reshape(pre_model.edge_pot,[],1);
    epot_all{n,2} = reshape(post_model.edge_pot,[],1);
    epot_in{n,1} = reshape(pre_model.edge_pot(Stim_cells,Stim_cells),[],1);
    epot_in{n,2} = reshape(post_model.edge_pot(Stim_cells,Stim_cells),[],1);
    epot_all_cum{n,1} = calc_cum_dist(epot_all{n,1},epot_bin_range);
    epot_all_cum{n,2} = calc_cum_dist(epot_all{n,2},epot_bin_range);
    epot_in_cum{n,1} = calc_cum_dist(epot_in{n,1},epot_bin_range);
    epot_in_cum{n,2} = calc_cum_dist(epot_in{n,2},epot_bin_range);
    
    % density
    dens_all{n,1} = sum(pre_model.graph(:))/num_node/(num_node-1);
    dens_all{n,2} = sum(post_model.graph(:))/num_node/(num_node-1);
    dens_in{n,1} = sum(sum(pre_model.graph(Stim_cells,Stim_cells)))/...
        num_stim_cell/(num_stim_cell-1);
    dens_in{n,2} = sum(sum(post_model.graph(Stim_cells,Stim_cells)))/...
        num_stim_cell/(num_stim_cell-1);
    
    % node degree
    ndeg_all{n,1} = full(sum(pre_model.graph,2)/2/num_node);
    ndeg_all{n,2} = full(sum(post_model.graph,2)/2/num_node);
    ndeg_in{n,1} = full(sum(pre_model.graph(Stim_cells,Stim_cells),2)/2/num_stim_cell);
    ndeg_in{n,2} = full(sum(post_model.graph(Stim_cells,Stim_cells),2)/2/num_stim_cell);
    ndeg_all_cum{n,1} = calc_cum_dist(ndeg_all{n,1},ndeg_bin_range);
    ndeg_all_cum{n,2} = calc_cum_dist(ndeg_all{n,2},ndeg_bin_range);
    ndeg_in_cum{n,1} = calc_cum_dist(ndeg_in{n,1},ndeg_bin_range);
    ndeg_in_cum{n,2} = calc_cum_dist(ndeg_in{n,2},ndeg_bin_range);
    
    % lcc
    lcc_all{n,1} = local_cluster_coeff(pre_model.graph);
    lcc_all{n,2} = local_cluster_coeff(post_model.graph);
    lcc_in{n,1} = local_cluster_coeff(pre_model.graph(Stim_cells,Stim_cells));
    lcc_in{n,2} = local_cluster_coeff(post_model.graph(Stim_cells,Stim_cells));
    lcc_all_cum{n,1} = calc_cum_dist(lcc_all{n,1},lcc_bin_range);
    lcc_all_cum{n,2} = calc_cum_dist(lcc_all{n,2},lcc_bin_range);
    lcc_in_cum{n,1} = calc_cum_dist(lcc_in{n,1},lcc_bin_range);
    lcc_in_cum{n,2} = calc_cum_dist(lcc_in{n,2},lcc_bin_range);
    
    % centrality
    cent_all{n,1} = eigenvec_centrality(pre_model.graph);
    cent_all{n,2} = eigenvec_centrality(post_model.graph);
    cent_in{n,1} = eigenvec_centrality(pre_model.graph(Stim_cells,Stim_cells));
    cent_in{n,2} = eigenvec_centrality(post_model.graph(Stim_cells,Stim_cells));
    cent_all_cum{n,1} = calc_cum_dist(cent_all{n,1},cent_bin_range);
    cent_all_cum{n,2} = calc_cum_dist(cent_all{n,2},cent_bin_range);
    cent_in_cum{n,1} = calc_cum_dist(cent_in{n,1},cent_bin_range);
    cent_in_cum{n,2} = calc_cum_dist(cent_in{n,2},cent_bin_range);
    
    % maximal cliques
    mc = maximalCliques(pre_model.graph);
    mc_pre = cell(size(mc,2),1);
    for ii = 1:size(mc,2)
        mc_pre{ii} = find(mc(:,ii));
    end
    mc = maximalCliques(post_model.graph);
    mc_post = cell(size(mc,2),1);
    for ii = 1:size(mc,2)
        mc_post{ii} = find(mc(:,ii));
    end
    
    % mc number
    mc_num(n,1) = length(mc_pre);
    mc_num(n,2) = length(mc_post);
    
    % mc size
    mc_sz{n,1} = cellfun('length',mc_pre);
    mc_sz{n,2} = cellfun('length',mc_post);
    mc_sz_cum{n,1} = calc_cum_dist(mc_sz{n,1},mc_sz_bin_range);
    mc_sz_cum{n,2} = calc_cum_dist(mc_sz{n,2},mc_sz_bin_range);
    
    %% core cells
    % cell with top 30% node degree
    core_ndeg = find(ndeg_all{n,2}>quantile(ndeg_all{n,2},0.8));
    
    % core with mc overlap
    mc_ov = cellfun(@(x) length(intersect(x,Stim_cells)),mc_post);
    keep_indx = mc_ov>=mc_sz{n,2}-1 & mc_sz{n,2}>=4;
    core_mc = unique(cell2mat(mc_post(keep_indx)));
    
    figure;set(gcf,'color','w');
    subplot(1,2,1)
    plotGraphHighlight(Coord_active,core_ndeg,mycc.orange,1);
    title('ndeg core')
    subplot(1,2,2)
    plotGraphHighlight(Coord_active,core_mc,mycc.orange,1);
    title('mc core')
    
end

save([save_path 'opto_spont_prop.mat'],'-v7.3');

%% plot node potentials and edge potentials
stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

figure; set(gcf,'color','w','position',[1984 320 430 463])

% mean node potentials
subplot(2,2,1); hold on
h = boxplot(cellfun(@(x) mean(x), npot_all(:,1)),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
h = boxplot(cellfun(@(x) mean(x), npot_all(:,2)),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([min(cell2mat(npot_all(:))),max(cell2mat(npot_all(:)))])
set(gca,'xtick',stepsz*[1,2],'xticklabel',{'pre','post'})
ylabel('node potential')

subplot(2,2,2); hold on
h = boxplot(cellfun(@(x) mean(x), npot_in(:,1)),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
h = boxplot(cellfun(@(x) mean(x), npot_in(:,2)),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([min(cell2mat(npot_in(:))),max(cell2mat(npot_in(:)))])
set(gca,'xtick',stepsz*[1,2],'xticklabel',{'pre','post'})
ylabel('node potential')

% mean edge potentials
subplot(2,2,3); hold on
h = boxplot(cellfun(@(x) mean(x), epot_all(:,1)),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
h = boxplot(cellfun(@(x) mean(x), epot_all(:,2)),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([min(cell2mat(epot_all(:))),max(cell2mat(epot_all(:)))])
set(gca,'xtick',stepsz*[1,2],'xticklabel',{'pre','post'})
ylabel('edge potential')

subplot(2,2,4); hold on
h = boxplot(cellfun(@(x) mean(x), epot_in(:,1)),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
h = boxplot(cellfun(@(x) mean(x), epot_in(:,2)),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([min(cell2mat(epot_in(:))),max(cell2mat(epot_in(:)))])
set(gca,'xtick',stepsz*[1,2],'xticklabel',{'pre','post'})
ylabel('edge potential')

saveas(gcf,[fig_path 'opto_spont_ep_np_prop.pdf']);

%% graph properties - whole model
figure; set(gcf,'color','w','position',[2454,301,560,465])

% density
boxwd = 0.2;
subplot(2,3,1);hold on;
dens_pre = cell2mat(dens_all(:,1)');
h = boxplot(dens_pre,'positions',0.5,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
dens_post = cell2mat(dens_all(:,2)');
h = boxplot(dens_post,'positions',1,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 1.5])
ylim([min([dens_pre,dens_post])-0.02 max([dens_pre,dens_post])]+0.02)
gcapos = get(gca,'position');
title('density')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);

% node degree
subplot(2,3,2)
plot_opto_spont_cum_hist(ndeg_all_cum,mycc,ndeg_bin_range,linew);
gcapos = get(gca,'position');
title('node degree');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% lcc
subplot(2,3,3)
plot_opto_spont_cum_hist(lcc_all_cum,mycc,lcc_bin_range,linew);
gcapos = get(gca,'position');
title('lcc');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% centrality
subplot(2,3,4)
plot_opto_spont_cum_hist(cent_all_cum,mycc,cent_bin_range,linew);
gcapos = get(gca,'position');
title('centrality');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% maximal clique number
boxwd = 0.2;
subplot(2,3,5);hold on;
h = boxplot(mc_num(:,1),'positions',0.5,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
h = boxplot(mc_num(:,2),'positions',1,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 1.5])
ylim([min(mc_num(:))-50 max(mc_num(:))+50])
gcapos = get(gca,'position');
title('NMC')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);

% maximal clique size
subplot(2,3,6)
plot_opto_spont_cum_hist(mc_sz_cum,mycc,mc_sz_bin_range,linew);
gcapos = get(gca,'position');
title('sMC');ylabel('p');
set(gca,'position',gcapos);
box on

saveas(gcf,[fig_path 'opto_spont_graph_prop_all.pdf']);

%% graph properties - stim network
figure; set(gcf,'color','w','position',[2183,140,377,469])

% density
boxwd = 0.2;
subplot(2,2,1);hold on;
dens_pre = full(cell2mat(dens_in(:,1)'));
h = boxplot(dens_pre,'positions',0.5,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
dens_post = full(cell2mat(dens_in(:,2)'));
h = boxplot(dens_post,'positions',1,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 1.5])
ylim([min([dens_pre;dens_post])-0.02 max([dens_pre;dens_post])]+0.02)
gcapos = get(gca,'position');
title('density')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);

% node degree
subplot(2,2,2)
plot_opto_spont_cum_hist(ndeg_in_cum,mycc,ndeg_bin_range,linew);
gcapos = get(gca,'position');
title('node degree');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% lcc
subplot(2,2,3)
plot_opto_spont_cum_hist(lcc_in_cum,mycc,lcc_bin_range,linew);
gcapos = get(gca,'position');
title('lcc');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% centrality
subplot(2,2,4)
plot_opto_spont_cum_hist(cent_in_cum,mycc,cent_bin_range,linew);
gcapos = get(gca,'position');
title('centrality');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

saveas(gcf,[fig_path 'opto_spont_graph_prop_in.pdf']);

end