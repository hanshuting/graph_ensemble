function [] = fig4_cc_crf_mc_graph_prop(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
k = param.k;
p = param.p;
cc_type = param.cc_type;
comm_type = param.comm_type;
ge_meth = param.ge_meth;
fig_path = param.fig_path.graph_prop;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;

ndeg_bin_range = param.ndeg_bin_range;
lcc_bin_range = param.lcc_bin_range;
cent_bin_range = param.cent_bin_range;
mc_sz_bin_range = param.mc_sz_bin_range;

load(ccode_path);

% initialize results
mc_sz = struct();
mc_num = struct();
dens = struct();
lcc = struct();
ndeg = struct();
cent = struct();

mc_num_cum = struct();
lcc_cum = struct();
ndeg_cum = struct();
cent_cum = struct();

%% go over experiments
expt_count = 0;
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        expt_count = expt_count+1;
        fprintf('Processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set path
        model_path = [result_path_base expt_name{n} '\models\']; 
        cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\']; 
        mc_path = [result_path_base expt_name{n} '\mc\'];
        rand_path = [result_path_base expt_name{n} '\rand_graph_' ge_meth '\'];
        
        % load data
        best_model = load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([cc_path  expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        load([mc_path expt_name{n} '_' expt_ee{e} '_mc_crf.mat']); % mc_crf
        load([mc_path expt_name{n} '_' expt_ee{e} '_mc_cc.mat']); % mc_cc
        load([rand_path expt_name{n} '_' expt_ee{e} '_mc_cc_rand.mat']);
        load([rand_path expt_name{n} '_' expt_ee{e} '_mc_crf_rand.mat']); 
        
        crf_graph = best_model.graph;
        edge_pot = best_model.edge_pot;
        num_node = size(crf_graph,1);
        
        %% plot graphs with latent neuron
        figure; set(gcf,'color','w','position',[600 600 900 300],'paperpositionmode','manual');
        h = subplot(1,2,1);
        plotLatentNeuronGraph(cc_graph,Coord_active,h);
        title('CC')
        h = subplot(1,2,2);
        plotLatentNeuronGraph(crf_graph,Coord_active,h);
        title('CRF')
        
        %% mc number
        mc_num(expt_count).cc = length(mc_cc);
        mc_num(expt_count).crf = length(mc_crf);
        mc_num(expt_count).cc_rand = cellfun('length',mc_cc_rand);
        mc_num(expt_count).crf_rand = cellfun('length',mc_crf_rand);
        
        %% mc size
        mc_sz(expt_count).crf = cellfun('length',mc_crf);
        mc_sz(expt_count).cc = cellfun('length',mc_cc);
        mc_sz(expt_count).crf_rand = cellfun('length',vertcat(mc_crf_rand{:}));
        mc_sz(expt_count).cc_rand = cellfun('length',vertcat(mc_cc_rand{:}));
        
        % cumulative distribution
        mc_sz_cum(expt_count).crf = calc_cum_dist(mc_sz(expt_count).crf,mc_sz_bin_range);
        mc_sz_cum(expt_count).cc = calc_cum_dist(mc_sz(expt_count).cc,mc_sz_bin_range);
        mc_sz_cum(expt_count).crf_rand = calc_cum_dist(mc_sz(expt_count).crf_rand,mc_sz_bin_range);
        mc_sz_cum(expt_count).cc_rand = calc_cum_dist(mc_sz(expt_count).cc_rand,mc_sz_bin_range);
        
        %% ------------- graph properties ---------------- %
        % 1. connection density
        dens(expt_count).cc = sum(cc_graph(:))/num_node/(num_node-1);
        dens(expt_count).crf = sum(crf_graph(:))/num_node/(num_node-1);
        dens(expt_count).cc_rand = sum(cc_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        dens(expt_count).crf_rand = sum(crf_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        
        % 2. node degree
        ndeg(expt_count).cc = sum(cc_graph,2)/2/num_node;
        ndeg(expt_count).crf = sum(crf_graph,2)/2/num_node;
        ndeg(expt_count).cc_rand = zeros(num_node,num_shuff);
        ndeg(expt_count).crf_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            ndeg(expt_count).cc_rand(:,i) = sum(squeeze(cc_rand_graph(:,:,i)),2)/2/num_node;
            ndeg(expt_count).crf_rand(:,i) = sum(squeeze(crf_rand_graph(:,:,i)),2)/2/num_node;
        end
        
        % distribution
        ndeg_cum(expt_count).cc = calc_cum_dist(ndeg(expt_count).cc,ndeg_bin_range);
        ndeg_cum(expt_count).crf = calc_cum_dist(ndeg(expt_count).crf,ndeg_bin_range);
        ndeg_cum(expt_count).cc_rand = calc_cum_dist(ndeg(expt_count).cc_rand,ndeg_bin_range);
        ndeg_cum(expt_count).crf_rand = calc_cum_dist(ndeg(expt_count).crf_rand,ndeg_bin_range);

        % 3. local clustering coefficient
        lcc(expt_count).cc = local_cluster_coeff(cc_graph);
        lcc(expt_count).crf = local_cluster_coeff(crf_graph);
        lcc(expt_count).cc_rand = zeros(num_node,num_shuff);
        lcc(expt_count).crf_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            lcc(expt_count).cc_rand(:,i) = local_cluster_coeff(cc_rand_graph(:,:,i));
            lcc(expt_count).crf_rand(:,i) = local_cluster_coeff(crf_rand_graph(:,:,i));
        end
        
        lcc_cum(expt_count).cc = calc_cum_dist(lcc(expt_count).cc,lcc_bin_range);
        lcc_cum(expt_count).crf = calc_cum_dist(lcc(expt_count).crf,lcc_bin_range);
        lcc_cum(expt_count).cc_rand = calc_cum_dist(lcc(expt_count).cc_rand,lcc_bin_range);
        lcc_cum(expt_count).crf_rand = calc_cum_dist(lcc(expt_count).crf_rand,lcc_bin_range);

        % 4. centrality
        cent(expt_count).crf = eigenvec_centrality(crf_graph);
        cent(expt_count).cc = eigenvec_centrality(cc_graph);
        cent(expt_count).crf_rand = zeros(num_node,num_shuff);
        cent(expt_count).cc_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            cent(expt_count).crf_rand(:,i) = eigenvec_centrality(crf_rand_graph(:,:,i));
            cent(expt_count).cc_rand(:,i) = eigenvec_centrality(cc_rand_graph(:,:,i));
        end
        
        cent_cum(expt_count).crf = calc_cum_dist(cent(expt_count).crf,cent_bin_range);
        cent_cum(expt_count).cc = calc_cum_dist(cent(expt_count).cc,cent_bin_range);
        cent_cum(expt_count).crf_rand = calc_cum_dist(cent(expt_count).crf_rand,cent_bin_range);
        cent_cum(expt_count).cc_rand = calc_cum_dist(cent(expt_count).cc_rand,cent_bin_range);
        
    end

end

%% save results
save([save_path 'graph_prop_mc_' savestr '.mat'],'expt_name','ee','mc_sz','mc_sz_cum',...
    'mc_num','dens','ndeg','ndeg_cum','lcc','lcc_cum','-v7.3');

%%
linew = 1;

figure;
set(gcf,'color','w','position',[1987 295 518 725],'PaperPositionMode','auto');

% density
boxwd = 0.2;
subplot(3,2,1);hold on;
h = boxplot(cell2mat(getNestedField(dens,'cc')),'positions',0.5,'width',...
    boxwd,'colors',mycc.green);
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
h = boxplot(cell2mat(getNestedField(dens,'crf')),'positions',1,'width',...
    boxwd,'colors',mycc.orange);
set(h(7,:),'visible','off')
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
xlim([0 1.5])
gcapos = get(gca,'position');
title('density')
set(gca,'xtick',[0.5 1],'xticklabel',{'CC','CRF'},'linewidth',linew)
set(gca,'position',gcapos);

% local clustering coefficient
subplot(3,2,2);
plot_graph_prop_cum_single(lcc_cum,mycc,lcc_bin_range);
gcapos = get(gca,'position');
title('clustering coeff');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% node degree
subplot(3,2,3);
plot_graph_prop_cum_single(ndeg_cum,mycc,ndeg_bin_range);
gcapos = get(gca,'position');
xlabel('node degree');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% centrality
subplot(3,2,4)
plot_graph_prop_cum_single(cent_cum,mycc,cent_bin_range);
gcapos = get(gca,'position');
title('centrality');ylabel('p');
set(gca,'position',gcapos);
legend off; box on

% mc number
boxwd = 0.2;
subplot(3,2,5);hold on;
h = boxplot(cell2mat(getNestedField(mc_num,'cc_rand')),'positions',0.5,'width',...
    boxwd,'colors',mycc.green_light);
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
h = boxplot(cell2mat(getNestedField(mc_num,'cc')),'positions',1,'width',...
    boxwd,'colors',mycc.green);
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
h = boxplot(cell2mat(getNestedField(mc_num,'crf_rand')),'positions',1.5,'width',...
    boxwd,'colors',mycc.orange_light);
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
h = boxplot(cell2mat(getNestedField(mc_num,'crf')),'positions',2,'width',...
    boxwd,'colors',mycc.orange);
set(h(7,:),'visible','off');set(h,'linewidth',linew);set(h(6,:),'linewidth',2*linew)
xlim([0 2.5])
gcapos = get(gca,'position');
title('NMC')
set(gca,'xtick',0.5:0.5:2,'xticklabel',{'CCrand','CC','CRFrand','CRF'},'linewidth',linew)
set(gca,'position',gcapos);

% mc size
subplot(3,2,6)
plot_graph_prop_cum_single(mc_sz_cum,mycc,mc_sz_bin_range);
gcapos = get(gca,'position');
title('sMC'); ylabel('p')
set(gca,'position',gcapos);
box on

% save figure
saveas(gcf,[fig_path 'graph_prop_mc_' savestr '.pdf']);

end
