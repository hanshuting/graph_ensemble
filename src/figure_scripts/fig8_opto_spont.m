function [] = fig8_opto_spont(param)

% parameters
expt_name = param.expt_name;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.opto_spont;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
linew = param.linew;
p = param.p;

ndeg_bin_range = param.ndeg_bin_range;

load(ccode_path);
graymap = load(param.graymap); % variable: cmap
bmap = load(param.bluemap); % variable: cmap

qnoise = 0.7;

%% initialize
num_expt = length(expt_name);

ndeg = cell(num_expt,2);
ndeg_in = cell(num_expt,2);
ndeg_out = cell(num_expt,2);

ndeg_dist = zeros(num_expt,length(ndeg_bin_range),2);
ndeg_in_dist = zeros(num_expt,length(ndeg_bin_range),2);
ndeg_out_dist = zeros(num_expt,length(ndeg_bin_range),2);

dens = cell(num_expt,2);
dens_in = cell(num_expt,2);
dens_out = cell(num_expt,2);

lcc = cell(num_expt,2);
lcc_in = cell(num_expt,2);
lcc_out = cell(num_expt,2);

cent = cell(num_expt,2);
cent_in = cell(num_expt,2);
cent_out = cell(num_expt,2);

ep_sum = cell(num_expt,1);
ep_sum_in = cell(num_expt,1);
ep_sum_out = cell(num_expt,1);

% number of edges for high-ranked neurons
pt_nedge = cell(num_expt,2);
pt_nedge_in = cell(num_expt,2);

%% collect data from all experiments
for n = 1:num_expt
    
    expt_ee = param.ee{n};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Stim_cells.mat']);
    
    pre_model = load([model_path expt_name{n} '_' expt_ee{1} ...
        '_loopy_best_model_' ge_type '.mat']);
    post_model = load([model_path expt_name{n} '_' expt_ee{2} ...
        '_loopy_best_model_' ge_type '.mat']);
    post_model.graph = full(post_model.graph);
    
    num_node = size(pre_model.graph,1);
    num_stim_cell = length(Stim_cells);
    nostim_cells = setdiff(1:num_node,Stim_cells);
    num_nostim_cell = num_node-num_stim_cell;
    num_neuron = num_node-1; % exclude added node
    
    pre_data = load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{1} '.mat']);
    pre_data = pre_data.data; %(:,1:num_neuron);
    post_data = load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{2} '.mat']);
    post_data = post_data.data; %(:,1:num_neuron);
    num_frame_pre = size(pre_data,1);
    num_frame_post = size(post_data,1);
    
    load([data_path expt_name{n} '\opto_stim_pre.mat']); % opto_stim_pre
    load([data_path expt_name{n} '\opto_stim_post.mat']); % opto_stim_post
    
    % convert to on edges
    pre_model.ep_on = getOnEdgePot(pre_model.graph,pre_model.G)';
    post_model.ep_on = getOnEdgePot(post_model.graph,post_model.G)';
    
    % shuffled models
    pre_shuffle = load([model_path 'shuffled_' expt_name{n} '_' ...
        expt_ee{1} '_loopy.mat']);
    for ii = 1:length(pre_shuffle.graphs)
        pre_shuffle.ep_on{ii} = getOnEdgePot(pre_shuffle.graphs{ii},...
            pre_shuffle.G{ii})';
        pre_shuffle.epsum{ii} = sum(pre_shuffle.ep_on{ii},2);
        pre_shuffle.epsum{ii}(sum(pre_shuffle.graphs{ii},2)==0) = NaN;
    end
    pre_shuffle.mepsum = nanmean(cellfun(@(x) nanmean(x),pre_shuffle.epsum));
    pre_shuffle.sdepsum = nanstd(cellfun(@(x) nanmean(x),pre_shuffle.epsum));
    post_shuffle = load([model_path 'shuffled_' expt_name{n} '_' ...
        expt_ee{2} '_loopy.mat']);
    for ii = 1:length(post_shuffle.graphs)
        post_shuffle.ep_on{ii} = getOnEdgePot(post_shuffle.graphs{ii},...
            post_shuffle.G{ii})';
        post_shuffle.epsum{ii} = sum(post_shuffle.ep_on{ii},2);
        post_shuffle.epsum{ii}(sum(post_shuffle.graphs{ii},2)==0) = NaN;
    end
    post_shuffle.mepsum = nanmean(cellfun(@(x) nanmean(x),post_shuffle.epsum));
    post_shuffle.sdepsum = nanstd(cellfun(@(x) nanmean(x),post_shuffle.epsum));
    
    % extend coordinates for add neuron model
    coords = Coord_active;
    coords(end+1,:) = [0 max(coords(:,2))];
    coords(end+1,:) = [0 0];
    
    %% visualize graph with circle layout, highlight pattern completion neuron
%     [~,hr_indx] = max(sum(post_model.ep_on(Stim_cells,Stim_cells),2));
%     hr_indx = Stim_cells(hr_indx); % 22 
    hr_indx = 22;
    nonpt_cell = setdiff(Stim_cells,hr_indx);
    num_half = round(num_nostim_cell/2); num_half(2) = num_nostim_cell-num_half;
    sort_indx = [nostim_cells(1:num_half(1)),hr_indx,...
        nostim_cells(num_half(1)+1:sum(num_half)),nonpt_cell];
    nodec = zeros(num_node,3);
    nodec(1:num_half(1),:) = repmat(mycc.gray_light,num_half(1),1);
    nodec(num_half(1)+1,:) = mycc.blue;
    nodec(num_half(1)+2:num_nostim_cell+1,:) = repmat(mycc.gray_light,num_half(2),1);
    nodec(num_nostim_cell+2:end,:) = repmat(mycc.blue,num_stim_cell-1,1);
    
    % set pre-model edge colors
    pre_edge_list = zeros(sum(pre_model.graph(:))/2,2);
    [pre_edge_list(:,2),pre_edge_list(:,1)] = find(tril(pre_model.graph(sort_indx,sort_indx)));
    pre_edgec = num2cell(pre_edge_list);
    for ii = 1:size(pre_edge_list,1)
        if ismember(pre_edge_list(ii,1),[num_half(1)+1,num_nostim_cell+2:num_node]) && ...
                ismember(pre_edge_list(ii,2),[num_half(1)+1,num_nostim_cell+2:num_node])
            pre_edgec{ii,3} = mycc.blue;
        else
            pre_edgec{ii,3} = mycc.gray_light; % NaN
        end
    end
    
    % set post-model edge colors
    post_edge_list = zeros(sum(post_model.graph(:))/2,2);
    [post_edge_list(:,2),post_edge_list(:,1)] = find(tril(post_model.graph(sort_indx,sort_indx)));
    post_edgec = num2cell(post_edge_list);
    for ii = 1:size(post_edge_list,1)
        if ismember(post_edge_list(ii,1),[num_half(1)+1,num_nostim_cell+2:num_node]) && ...
                ismember(post_edge_list(ii,2),[num_half(1)+1,num_nostim_cell+2:num_node])
            post_edgec{ii,3} = mycc.blue;
        else
            post_edgec{ii,3} = mycc.gray_light; % NaN;
        end
    end
    
    % plot
    figure; set(gcf,'color','w','position',[2022 313 894 372]);
    subplot(1,2,1); title('pre')
    visGraphCirc(pre_model.graph(sort_indx,sort_indx),'edgeColor',pre_edgec,...
        'nodeColor',nodec);
    subplot(1,2,2); title('post')
    visGraphCirc(post_model.graph(sort_indx,sort_indx),'edgeColor',post_edgec,...
        'nodeColor',nodec);
    
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_circlelayout_' ...
        ge_type '_stim_neuron.pdf']);
    
    %% plot pre and post models
    figure; set(gcf,'color','w','position',[2154 340 941 597])
    cc_range = [min(pre_model.ep_on(:)) max(pre_model.ep_on(:))];
    ep_range = [-1 0.2];
    subplot(1,2,1)
    plotGraphModelHighlightEP(pre_model.graph,coords,pre_model.ep_on,...
        cc_range,ep_range,graymap.cmap,[]);
    subplot(1,2,2)
%     plotGraphModelHighlightEP(post_model.graph,coords,post_model.ep_on,...
%         cc_range,ep_range,bmap.cmap,[]);
    plotGraphModelHighlightEP(post_model.graph,coords,post_model.ep_on,...
        cc_range,ep_range,graymap.cmap,[]);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_ON_unnormalized_' ...
        ge_type '_pre_post_models.pdf']);
    
    %% graph properties
    % edge potential sum
    ep_sum{n,1} = sum(pre_model.ep_on,2);
    ep_sum{n,2} = sum(post_model.ep_on,2);
    ep_sum_in{n,1} = ep_sum{n,1}(Stim_cells);
    ep_sum_in{n,2} = ep_sum{n,2}(Stim_cells);
    ep_sum_out{n,1} = ep_sum{n,1}(nostim_cells);
    ep_sum_out{n,2} = ep_sum{n,2}(nostim_cells);
    
    % density
    dens{n,1} = sum(sum(pre_model.graph))/num_node/(num_node-1);
    dens{n,2} = sum(sum(post_model.graph))/num_node/(num_node-1);
    dens_in{n,1} = sum(sum(pre_model.graph(Stim_cells,Stim_cells)))/...
        num_stim_cell/(num_stim_cell-1);
    dens_in{n,2} = sum(sum(post_model.graph(Stim_cells,Stim_cells)))/...
        num_stim_cell/(num_stim_cell-1);
    dens_out{n,1} = sum(sum(pre_model.graph(nostim_cells,nostim_cells)))/...
        num_nostim_cell/(num_nostim_cell-1);
    dens_out{n,2} = sum(sum(post_model.graph(nostim_cells,nostim_cells)))/...
        num_nostim_cell/(num_nostim_cell-1);
    
    % node degree
    ndeg{n,1} = sum(pre_model.graph,2)/2/num_node;
    ndeg{n,2} = sum(post_model.graph,2)/2/num_node;
    ndeg_in{n,1} = ndeg{n,1}(Stim_cells);
    ndeg_in{n,2} = ndeg{n,2}(Stim_cells);
    ndeg_out{n,1} = ndeg{n,1}(nostim_cells);
    ndeg_out{n,2} = ndeg{n,2}(nostim_cells);
    
    % node degree distribution
    ndeg_dist(n,:,1) = histc(ndeg{n,1},ndeg_bin_range);
    ndeg_dist(n,:,1) = ndeg_dist(n,:,1)/sum(ndeg_dist(n,:,1));
    ndeg_dist(n,:,2) = histc(ndeg{n,2},ndeg_bin_range);
    ndeg_dist(n,:,2) = ndeg_dist(n,:,2)/sum(ndeg_dist(n,:,2));
    ndeg_in_dist(n,:,1) = histc(ndeg_in{n,1},ndeg_bin_range);
    ndeg_in_dist(n,:,1) = ndeg_in_dist(n,:,1)/sum(ndeg_in_dist(n,:,1));
    ndeg_in_dist(n,:,2) = histc(ndeg_in{n,2},ndeg_bin_range);
    ndeg_in_dist(n,:,2) = ndeg_in_dist(n,:,2)/sum(ndeg_in_dist(n,:,2));
    ndeg_out_dist(n,:,1) = histc(ndeg_out{n,1},ndeg_bin_range);
    ndeg_out_dist(n,:,1) = ndeg_out_dist(n,:,1)/sum(ndeg_out_dist(n,:,1));
    ndeg_out_dist(n,:,2) = histc(ndeg_out{n,2},ndeg_bin_range);
    ndeg_out_dist(n,:,2) = ndeg_out_dist(n,:,2)/sum(ndeg_out_dist(n,:,2));
    
    % lcc
    lcc{n,1} = local_cluster_coeff(pre_model.graph);
    lcc{n,2} = local_cluster_coeff(post_model.graph);
    lcc_in{n,1} = lcc{n,1}(Stim_cells);
    lcc_in{n,2} = lcc{n,2}(Stim_cells);
    lcc_out{n,1} = lcc{n,1}(nostim_cells);
    lcc_out{n,2} = lcc{n,2}(nostim_cells);
    
    % centrality
    cent{n,1} = eigenvec_centrality(pre_model.graph);
    cent{n,2} = eigenvec_centrality(post_model.graph);
    cent_in{n,1} = cent{n,1}(Stim_cells);
    cent_in{n,2} = cent{n,2}(Stim_cells);
    cent_out{n,1} = cent{n,1}(nostim_cells);
    cent_out{n,2} = cent{n,2}(nostim_cells);
    
    %% predict with cosine similarity
    label_pre = opto_stim_pre'==1;
    label_post = opto_stim_post'==1;
    
    figure; set(gcf,'color','w','position',[2014 162 476 432])
    
    % pre
    subplot(2,2,1); hold on
    % random ensemble performance
    rd_ens = zeros(100,length(Stim_cells));
    auc_rd_pre = zeros(100,1);
    for ii = 1:100
        rd_ens(ii,:) = randperm(num_neuron,length(Stim_cells));
        [~,sim_core] = core_cos_sim(rd_ens(ii,:),pre_data(:,1:num_neuron),label_pre);
        [xx,yy,~,auc_rd_pre(ii)] = perfcurve(label_pre,sim_core,1);
    end
    plot(xx,yy,'color','k','linewidth',2*linew);
    % individual cells
    cc = jet(64);
    auc_pre = zeros(num_stim_cell,1);
    for ii = 1:num_stim_cell
        [~,sim_core] = core_cos_sim(Stim_cells(ii),pre_data(:,1:num_neuron),label_pre);
        [xx,yy,~,auc_pre(ii)] = perfcurve(label_pre,sim_core,1);
        curve_cc = cc(ceil((ep_sum{n,1}(Stim_cells(ii))-min(ep_sum{n,1}(Stim_cells)))/...
            (max(ep_sum{n,1}(Stim_cells))-min(ep_sum{n,1}(Stim_cells)))*63+1),:);
        plot(xx,yy,'color',curve_cc,'linewidth',linew);
    end
    xlabel('FPR'); ylabel('TPR'); title('pre')
    % plot auc
    nodesz = 30;
    subplot(2,2,3); hold on
    scatter(ep_sum{n,1}(Stim_cells),auc_pre,nodesz,mycc.gray,'filled')
    scatter(ep_sum{n,1}(hr_indx),auc_pre(Stim_cells==hr_indx),nodesz,mycc.red,'filled')
    xlabel('node strength'); ylabel('AUC');
    
    % post
    subplot(2,2,2); hold on
    % random ensemble performance
    auc_rd_post = zeros(100,1);
    for ii = 1:100
        [~,sim_core] = core_cos_sim(rd_ens(ii,:),post_data(:,1:num_neuron),label_post);
        [xx,yy,~,auc_rd_post(ii)] = perfcurve(label_post,sim_core,1);
    end
    plot(xx,yy,'color','k','linewidth',2*linew);
    % individual cells
    cc = jet(64);
    auc_post = zeros(num_stim_cell,1);
    for ii = 1:num_stim_cell
        [~,sim_core] = core_cos_sim(Stim_cells(ii),post_data(:,1:num_neuron),label_post);
        [xx,yy,~,auc_post(ii)] = perfcurve(label_post,sim_core,1);
        curve_cc = cc(ceil((ep_sum{n,2}(Stim_cells(ii))-min(ep_sum{n,2}(Stim_cells)))/...
            (max(ep_sum{n,2}(Stim_cells))-min(ep_sum{n,2}(Stim_cells)))*63+1),:);
        plot(xx,yy,'color',curve_cc,'linewidth',linew);
    end
    xlabel('FPR'); ylabel('TPR'); title('post')
    % plot auc
    nodesz = 30;
    subplot(2,2,4); hold on
    scatter(ep_sum{n,2}(Stim_cells),auc_post,nodesz,mycc.gray,'filled')
    scatter(ep_sum{n,2}(hr_indx),auc_post(Stim_cells==hr_indx),nodesz,mycc.red,'filled')
    xlabel('node strength'); ylabel('AUC');
    
    % plot ensemble auc, scale figures
    aucmi = min([auc_pre;auc_post])-0.1;
    aucma = max([auc_pre;auc_post])+0.1;
    nsmi = min([ep_sum{n,1}(Stim_cells)',ep_sum{n,2}(Stim_cells)']);
    nsma = max([ep_sum{n,1}(Stim_cells)',ep_sum{n,2}(Stim_cells)']);
    subplot(2,2,3); hold on
    plot([nsmi nsma],mean(auc_rd_pre)*[1 1],'k--')
    plot([nsmi nsma],(mean(auc_rd_pre)+std(auc_rd_pre))*[1 1],'--',...
        'color',mycc.gray_light);
    plot([nsmi nsma],(mean(auc_rd_pre)-std(auc_rd_pre))*[1 1],'--',...
        'color',mycc.gray_light);
    plot(pre_shuffle.mepsum*[1 1],[aucmi aucma],'k--');
    plot((pre_shuffle.mepsum+pre_shuffle.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',mycc.gray_light);
    plot((pre_shuffle.mepsum-pre_shuffle.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',mycc.gray_light);
    xlim([nsmi nsma]); ylim([aucmi aucma]);
    subplot(2,2,4); hold on
    plot([nsmi nsma],mean(auc_rd_post)*[1 1],'k--')
    plot([nsmi nsma],(mean(auc_rd_post)+std(auc_rd_post))*[1 1],'--',...
        'color',mycc.gray_light);
    plot([nsmi nsma],(mean(auc_rd_post)-std(auc_rd_post))*[1 1],'--',...
        'color',mycc.gray_light);
    plot(post_shuffle.mepsum*[1 1],[aucmi aucma],'k--');
    plot((post_shuffle.mepsum+post_shuffle.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',mycc.gray_light);
    plot((post_shuffle.mepsum-post_shuffle.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',mycc.gray_light);
    xlim([nsmi nsma]); ylim([aucmi aucma]);
    
    saveas(gcf,[fig_path 'opto_spont_ROC.pdf']);
    
    % cells that can do PT (index in Stim_cells)
    pt_indx = find(auc_post>=(mean(auc_rd_post)+std(auc_rd_post)) & ...
        ep_sum{n,2}(Stim_cells)>=(post_shuffle.mepsum+post_shuffle.sdepsum));
    
    % high-ranked neuron connections
    pt_nedge{n,1} = sum(pre_model.graph(Stim_cells(pt_indx),:),2);
    pt_nedge{n,2} = sum(post_model.graph(Stim_cells(pt_indx),:),2);
    pt_nedge_in{n,1} = sum(pre_model.graph(Stim_cells(pt_indx),Stim_cells),2);
    pt_nedge_in{n,2} = sum(post_model.graph(Stim_cells(pt_indx),Stim_cells),2);
    
    
    %% predict with model
    % change single node activity and predict with LL
%     LL_frame_pre_on = zeros(num_stim_cell,num_frame_pre,2);
%     LL_frame_pre_off = zeros(num_stim_cell,num_frame_pre,2);
%     LL_frame_post_on = zeros(num_stim_cell,num_frame_post,2);
%     LL_frame_post_off = zeros(num_stim_cell,num_frame_post,2);
%     for ii = 1:num_stim_cell
%         for jj = 1:num_frame_pre
%             frame_vec = pre_data(jj,:);
%             % stim on
%             frame_vec(end) = 1;
%             frame_vec(Stim_cells(ii)) = 0;
%             LL_frame_pre_on(ii,jj,1) = compute_avg_log_likelihood(pre_model.node_pot,...
%                 pre_model.edge_pot,pre_model.logZ,frame_vec);
%             frame_vec(Stim_cells(ii)) = 1;
%             LL_frame_pre_on(ii,jj,2) = compute_avg_log_likelihood(pre_model.node_pot,...
%                 pre_model.edge_pot,pre_model.logZ,frame_vec);
%             % stim off
%             frame_vec(end) = 0;
%             frame_vec(Stim_cells(ii)) = 0;
%             LL_frame_pre_off(ii,jj,1) = compute_avg_log_likelihood(pre_model.node_pot,...
%                 pre_model.edge_pot,pre_model.logZ,frame_vec);
%             frame_vec(Stim_cells(ii)) = 1;
%             LL_frame_pre_off(ii,jj,2) = compute_avg_log_likelihood(pre_model.node_pot,...
%                 pre_model.edge_pot,pre_model.logZ,frame_vec);
%         end
%         for jj = 1:num_frame_post
%             frame_vec = post_data(jj,:);
%             % stim on
%             frame_vec(end) = 1;
%             frame_vec(Stim_cells(ii)) = 0;
%             LL_frame_post_on(ii,jj,1) = compute_avg_log_likelihood(post_model.node_pot,...
%                 post_model.edge_pot,post_model.logZ,frame_vec);
%             frame_vec(Stim_cells(ii)) = 1;
%             LL_frame_post_on(ii,jj,2) = compute_avg_log_likelihood(post_model.node_pot,...
%                 post_model.edge_pot,post_model.logZ,frame_vec);
%             % stim off
%             frame_vec(end) = 0;
%             frame_vec(Stim_cells(ii)) = 0;
%             LL_frame_post_off(ii,jj,1) = compute_avg_log_likelihood(post_model.node_pot,...
%                 post_model.edge_pot,post_model.logZ,frame_vec);
%             frame_vec(Stim_cells(ii)) = 1;
%             LL_frame_post_off(ii,jj,2) = compute_avg_log_likelihood(post_model.node_pot,...
%                 post_model.edge_pot,post_model.logZ,frame_vec);
%         end
%     end
%     LL_stim_pre = squeeze(LL_frame_pre_on(:,:,2)-LL_frame_pre_on(:,:,1));
%     LL_nostim_pre = squeeze(LL_frame_pre_off(:,:,2)-LL_frame_pre_off(:,:,1));
%     LL_rel_pre = LL_stim_pre-LL_nostim_pre;
%     LL_stim_post = squeeze(LL_frame_post_on(:,:,2)-LL_frame_post_on(:,:,1));
%     LL_nostim_post = squeeze(LL_frame_post_off(:,:,2)-LL_frame_post_off(:,:,1));
%     LL_rel_post = LL_stim_post-LL_nostim_post;
%     
%     % plot
%     figure; set(gcf,'color','w','position',[1979 430 485 430])
%     cc = jet(64); nodesz = 30;
%     % ROC - pre
%     subplot(2,2,1); hold on
%     auc_pre = zeros(num_stim_cell,1);
%     for ii = 1:num_stim_cell
%         [xx,yy,~,auc_pre(ii)] = perfcurve(opto_stim_pre'==1,LL_stim_pre(ii,:),1);
%         curve_cc = cc(ceil((ep_sum{n,1}(Stim_cells(ii))-min(ep_sum{n,1}(Stim_cells)))/...
%             (max(ep_sum{n,1}(Stim_cells))-min(ep_sum{n,1}(Stim_cells)))*63+1),:);
%         plot(xx,yy,'color',curve_cc);
%     end
%     xlabel('FPR'); ylabel('TPR'); title('pre')
%     % ROC - post
%     subplot(2,2,2); hold on
%     auc_post = zeros(num_stim_cell,1);
%     for ii = 1:num_stim_cell
%         [xx,yy,~,auc_post(ii)] = perfcurve(opto_stim_post'==1,LL_stim_post(ii,:),1);
%         curve_cc = cc(ceil((ep_sum{n,2}(Stim_cells(ii))-min(ep_sum{n,2}(Stim_cells)))/...
%             (max(ep_sum{n,2}(Stim_cells))-min(ep_sum{n,2}(Stim_cells)))*63+1),:);
%         plot(xx,yy,'color',curve_cc);
%     end
%     xlabel('FPR'); ylabel('TPR'); title('post')
%     % AUC vs NS - pre
%     subplot(2,2,3); hold on
%     scatter(ep_sum{n,1}(Stim_cells),auc_pre,nodesz,mycc.gray,'filled')
%     scatter(ep_sum{n,1}(hr_indx),auc_pre(Stim_cells==hr_indx),nodesz,mycc.red,'filled')
%     xlabel('node strength'); ylabel('AUC');
%     % AUC vs NS - post
%     subplot(2,2,4); hold on
%     scatter(ep_sum{n,2}(Stim_cells),auc_post,nodesz,mycc.gray,'filled')
%     scatter(ep_sum{n,2}(hr_indx),auc_post(Stim_cells==hr_indx),nodesz,mycc.red,'filled')
%     xlabel('node strength'); ylabel('AUC');
%     
%     
end

% save([save_path 'opto_spont_prop.mat'],'-v7.3');

%% distribution of P(k)
% nsz = 5;
% figure;
% subplot(2,3,1); hold on
% for ii = 1:num_expt
%     plot(ndeg_bin_range,ndeg_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
%     plot(ndeg_bin_range,ndeg_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
%     scatter(ndeg_bin_range,ndeg_dist(ii,:,1),nsz,mycc.gray,'o','filled')
%     scatter(ndeg_bin_range,ndeg_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
% end
% 
% subplot(2,3,2); hold on
% for ii = 1:num_expt
%     plot(ndeg_bin_range,ndeg_in_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
%     plot(ndeg_bin_range,ndeg_in_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
%     scatter(ndeg_bin_range,ndeg_in_dist(ii,:,1),nsz,mycc.gray,'o','filled')
%     scatter(ndeg_bin_range,ndeg_in_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
% end
% 
% subplot(2,3,3); hold on
% for ii = 1:num_expt
%     plot(ndeg_bin_range,ndeg_out_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
%     plot(ndeg_bin_range,ndeg_out_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
%     scatter(ndeg_bin_range,ndeg_out_dist(ii,:,1),nsz,mycc.gray,'o','filled')
%     scatter(ndeg_bin_range,ndeg_out_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
% end
% 
% subplot(2,3,4); hold on
% for ii = 1:num_expt
%     plot(log(ndeg_bin_range),log(ndeg_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
%     plot(log(ndeg_bin_range),log(ndeg_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
%     scatter(log(ndeg_bin_range),log(ndeg_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
%     scatter(log(ndeg_bin_range),log(ndeg_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
% end
% 
% subplot(2,3,5); hold on
% for ii = 1:num_expt
%     plot(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
%     plot(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
%     scatter(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
%     scatter(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
% end
% 
% subplot(2,3,6); hold on
% for ii = 1:num_expt
%     plot(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
%     plot(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
%     scatter(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
%     scatter(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
% end

%% number of connections of PT/nonPT cells
boxwd = 0.2;
stepsz = 0.5;

% stim network
figure; set(gcf,'color','w'); 
subplot(1,2,1); hold on
nedge_pre = cell2mat(pt_nedge_in(:,1)');
nedge_post = cell2mat(pt_nedge_in(:,2)');
h = boxplot(nedge_pre,'positions',stepsz,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
h = boxplot(nedge_post,'positions',2*stepsz,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 3*stepsz]);
ylim([min([nedge_pre;nedge_post])-0.02 max([nedge_pre;nedge_post])]+0.02)
gcapos = get(gca,'position');
ylabel('# connections')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);
[~,pval] = ttest2(nedge_pre,nedge_post);
title(num2str(pval));
box off

% whole network
subplot(1,2,2); hold on
nedge_pre = cell2mat(pt_nedge(:,1)');
nedge_post = cell2mat(pt_nedge(:,2)');
h = boxplot(nedge_pre,'positions',stepsz,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
h = boxplot(nedge_post,'positions',2*stepsz,'width',boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
setBoxStyle(h,linew);
xlim([0 3*stepsz]);
ylim([min([nedge_pre;nedge_post])-0.02 max([nedge_pre;nedge_post])]+0.02)
gcapos = get(gca,'position');
ylabel('# connections')
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);
[~,pval] = ttest2(nedge_pre,nedge_post);
title(num2str(pval));
box off

saveas(gcf,[fig_path 'pt_cell_connections.pdf']);

%% graph properties - whole network
stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

figure; set(gcf,'color','w','position',[2017 597 717 148])

% density
boxwd = 0.2;
subplot(1,5,1);hold on;
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

% high-ranked neuron connections
% subplot(1,6,2);hold on;
% pp = plot_pair_graph(cell2mat(hr_nedge(:,1)),cell2mat(hr_nedge(:,2)),...
%     mycc.black,mycc.blue,p);
% gcapos = get(gca,'position');
% ylabel('high-ranked neuron connections'); title(num2str(pp));
% set(gca,'position',gcapos);
% legend off; box off

% edge pot sum
subplot(1,5,2)
pp = plot_pair_graph(cell2mat(ep_sum(:,1)),cell2mat(ep_sum(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node strengh'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(1,5,3)
pp = plot_pair_graph(cell2mat(ndeg(:,1)),cell2mat(ndeg(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node degree');title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(1,5,4)
pp = plot_pair_graph(cell2mat(lcc(:,1)),cell2mat(lcc(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(1,5,5)
pp = plot_pair_graph(cell2mat(cent(:,1)),cell2mat(cent(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('all')
saveas(gcf,[fig_path 'opto_spont_graph_prop_whole_unnormalized_on.pdf']);

%% stim network
figure; set(gcf,'color','w','position',[2017 597 717 148])

% density
boxwd = 0.2;
subplot(1,5,1);hold on;
dens_pre = cell2mat(dens_in(:,1)')*100;
h = boxplot(dens_pre,'positions',stepsz,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
dens_post = cell2mat(dens_in(:,2)')*100;
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

% edge pot sum
subplot(1,5,2)
pp = plot_pair_graph(cell2mat(ep_sum_in(:,1)),cell2mat(ep_sum_in(:,2)),...
    mycc.black,mycc.blue,p);
plot([0.4*ones(length(pt_indx),1) 0.6*ones(length(pt_indx),1)]',...
    [ep_sum_in{1,1}(pt_indx) ep_sum_in{1,2}(pt_indx)]','r');
gcapos = get(gca,'position'); title(num2str(pp));
ylabel('node strengh');
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(1,5,3)
pp = plot_pair_graph(cell2mat(ndeg_in(:,1)),cell2mat(ndeg_in(:,2)),...
    mycc.black,mycc.blue,p);
plot([0.4*ones(length(pt_indx),1) 0.6*ones(length(pt_indx),1)]',...
    [ndeg_in{1,1}(pt_indx) ndeg_in{1,2}(pt_indx)]','r');
gcapos = get(gca,'position');
ylabel('node degree'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(1,5,4)
pp = plot_pair_graph(cell2mat(lcc_in(:,1)),cell2mat(lcc_in(:,2)),...
    mycc.black,mycc.blue,p);
plot([0.4*ones(length(pt_indx),1) 0.6*ones(length(pt_indx),1)]',...
    [lcc_in{1,1}(pt_indx) lcc_in{1,2}(pt_indx)]','r');
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(1,5,5)
pp = plot_pair_graph(cell2mat(cent_in(:,1)),cell2mat(cent_in(:,2)),...
    mycc.black,mycc.blue,p);
plot([0.4*ones(length(pt_indx),1) 0.6*ones(length(pt_indx),1)]',...
    [cent_in{1,1}(pt_indx) cent_in{1,2}(pt_indx)]','r');
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('stim')
saveas(gcf,[fig_path 'opto_spont_graph_prop_stim_network_unnormalized_on.pdf']);

%% nostim network
figure; set(gcf,'color','w','position',[2017 597 717 148])

% density
boxwd = 0.2;
subplot(1,5,1);hold on;
dens_pre = cell2mat(dens_out(:,1)')*100;
h = boxplot(dens_pre,'positions',stepsz,'width',boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
dens_post = cell2mat(dens_out(:,2)')*100;
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

% edge pot sum
subplot(1,5,2)
pp = plot_pair_graph(cell2mat(ep_sum_out(:,1)),cell2mat(ep_sum_out(:,2)),...
    mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node strengh'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(1,5,3)
pp = plot_pair_graph(cell2mat(ndeg_out(:,1)),cell2mat(ndeg_out(:,2)),...
    mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node degree'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(1,5,4)
pp = plot_pair_graph(cell2mat(lcc_out(:,1)),cell2mat(lcc_out(:,2)),...
    mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(1,5,5)
pp = plot_pair_graph(cell2mat(cent_out(:,1)),cell2mat(cent_out(:,2)),...
    mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('no stim')
saveas(gcf,[fig_path 'opto_spont_graph_prop_nostim_network_unnormalized_on.pdf']);


end