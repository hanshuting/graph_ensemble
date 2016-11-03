function [] = figS2_compare_add_neuron(param)

% parameters
expt_name = param.expt_name;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.graph_prop;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
epsum_bin_range = param.epsum_bin_range;
linew = param.linew;
p = param.p;

ndeg_bin_range = param.ndeg_bin_range;
lcc_bin_range = param.lcc_bin_range;
cent_bin_range = param.cent_bin_range;

load(ccode_path);
graymap = load(param.graymap); % variable: cmap
bmap = load(param.bluemap); % variable: cmap

%% initialize
ndeg = cell(length(expt_name),2);
ndeg_cum = cell(length(expt_name),2);

dens = cell(length(expt_name),2);

lcc = cell(length(expt_name),2);
lcc_cum = cell(length(expt_name),2);

cent = cell(length(expt_name),2);
cent_cum = cell(length(expt_name),2);

ep_sum = cell(length(expt_name),1);
ep_sum_cum = cell(length(expt_name),1);

%% collect data from all experiments
for n = 1:length(expt_name)
    
    expt_ee = param.ee{n};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    base_model = load([model_path expt_name{n} '_' expt_ee{1} '_loopy_best_model_' ge_type '.mat']);
    an_model = load([model_path expt_name{n} '_' expt_ee{2} '_loopy_best_model_' ge_type '.mat']);
    
    % convert to on edges
    base_model.ep_on = getOnEdgePot(base_model.graph,base_model.G)';
    an_model.ep_on = getOnEdgePot(an_model.graph,an_model.G)';
    
    % extend coordinates for add neuron model
    coords = Coord_active;
    coords(end+1,:) = [0 max(coords(:,2))];
    coords(end+1,:) = [0 0];
    
    % plot pre and post models
    figure; set(gcf,'color','w','position',[2154 340 941 597])
    cc_range = [min(base_model.ep_on(:)) max(base_model.ep_on(:))];
    ep_range = [-1 0.2];
    subplot(1,2,1)
    plotGraphModelHighlightEP(base_model.graph,Coord_active,base_model.ep_on,...
        cc_range,ep_range,graymap.cmap,[]);
    subplot(1,2,2)
    plotGraphModelHighlightEP(an_model.graph,coords,an_model.ep_on,...
        cc_range,ep_range,bmap.cmap,[]);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} ...
        '_base_add_neuron_model_graph.pdf']);

    % edge potential sum
    ep_sum{n,1} = sum(base_model.ep_on,2);
    ep_sum{n,2} = sum(an_model.ep_on(1:num_node,1:num_node),2);
    ep_sum_cum{n,1} = calc_cum_dist(ep_sum{n,1},epsum_bin_range);
    ep_sum_cum{n,2} = calc_cum_dist(ep_sum{n,2},epsum_bin_range);
    
    % density
    dens{n,1} = sum(sum(base_model.graph))/num_node/(num_node-1);
    dens{n,2} = sum(sum(an_model.graph))/num_node/(num_node-1);
    
    % node degree
    ndeg{n,1} = sum(base_model.graph,2)/2/num_node;
    ndeg{n,2} = sum(an_model.graph(1:num_node,1:num_node),2)/2/num_node;
    ndeg_cum{n,1} = calc_cum_dist(ndeg{n,1},ndeg_bin_range);
    ndeg_cum{n,2} = calc_cum_dist(ndeg{n,2},ndeg_bin_range);
    
    % lcc
    lcc{n,1} = local_cluster_coeff(base_model.graph);
    lcc{n,2} = local_cluster_coeff(an_model.graph(1:num_node,1:num_node));
    lcc_cum{n,1} = calc_cum_dist(lcc{n,1},lcc_bin_range);
    lcc_cum{n,2} = calc_cum_dist(lcc{n,2},lcc_bin_range);
    
    % centrality
    cent{n,1} = eigenvec_centrality(base_model.graph);
    cent{n,2} = eigenvec_centrality(an_model.graph(1:num_node,1:num_node));
    cent_cum{n,1} = calc_cum_dist(cent{n,1},cent_bin_range);
    cent_cum{n,2} = calc_cum_dist(cent{n,2},cent_bin_range);

end


%% graph properties - nostim network
stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

figure; set(gcf,'color','w','position',[2006 320 691 425])

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
set(gca,'xtick',[0.5 1],'xticklabel',{'base','hidden'},'linewidth',linew)
set(gca,'position',gcapos);
box off

% edge pot sum
subplot(2,3,2)
% plot_opto_spont_cum_hist(ep_sum_cum,mycc,epsum_bin_range,linew);
plot_pair_graph(cell2mat(ep_sum(:,1)),cell2mat(ep_sum(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
% xlabel('sum(edge pot)');ylabel('log(p)');
ylabel('rank')
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(2,3,3)
% plot_opto_spont_cum_hist(ndeg_cum,mycc,ndeg_bin_range,linew);
plot_pair_graph(cell2mat(ndeg(:,1)),cell2mat(ndeg(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
% xlabel('node degree');ylabel('log(p)');
ylabel('node degree')
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(2,3,4)
% plot_opto_spont_cum_hist(lcc_cum,mycc,lcc_bin_range,linew);
plot_pair_graph(cell2mat(lcc(:,1)),cell2mat(lcc(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
% xlabel('clustering coeff');ylabel('log(p)');
ylabel('clustering coeff')
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(2,3,5)
% plot_opto_spont_cum_hist(cent_cum,mycc,cent_bin_range,linew);
plot_pair_graph(cell2mat(cent(:,1)),cell2mat(cent(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
% xlabel('centrality');ylabel('log(p)');
ylabel('centrality')
set(gca,'position',gcapos);
legend off; box off

saveas(gcf,[fig_path 'compare_base_add_neuron_graph_prop.pdf']);


end