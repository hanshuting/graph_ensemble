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
hr_nedge = cell(num_expt,2);

%% collect data from all experiments
for n = 1:num_expt
    
    expt_ee = param.ee{n};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Stim_cells.mat']);
    
    pre_model = load([model_path expt_name{n} '_' expt_ee{1} '_loopy_best_model_' ge_type '.mat']);
    post_model = load([model_path expt_name{n} '_' expt_ee{2} '_loopy_best_model_' ge_type '.mat']);
    post_model.graph = full(post_model.graph);
    
    num_node = size(pre_model.graph,1);
    num_stim_cell = length(Stim_cells);
    nostim_cells = setdiff(1:num_node,Stim_cells);
    num_nostim_cell = num_node-num_stim_cell;
    num_neuron = num_node-1; % exclude added node
    
    pre_data = load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{1} '.mat']);
    pre_data = pre_data.data(:,1:num_neuron);
    post_data = load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{2} '.mat']);
    post_data = post_data.data(:,1:num_neuron);
    num_frame_pre = size(pre_data,1);
    num_frame_post = size(post_data,1);
    
    load([data_path expt_name{n} '\opto_stim_pre.mat']); % opto_stim_pre
    load([data_path expt_name{n} '\opto_stim_post.mat']); % opto_stim_post
    
    % convert to on edges
    pre_model.ep_on = getOnEdgePot(pre_model.graph,pre_model.G)';
    post_model.ep_on = getOnEdgePot(post_model.graph,post_model.G)';
    
    % extend coordinates for add neuron model
    coords = Coord_active;
    coords(end+1,:) = [0 max(coords(:,2))];
    coords(end+1,:) = [0 0];
    
    %% visualize graph with circle layout, highlight pattern completion neuron
    [~,hr_indx] = max(sum(post_model.ep_on(Stim_cells,Stim_cells),2));
    hr_indx = Stim_cells(hr_indx); % 22 
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
    plotGraphModelHighlightEP(post_model.graph,coords,post_model.ep_on,...
        cc_range,ep_range,bmap.cmap,[]);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_ON_unnormalized_' ...
        ge_type '_pre_post_models.pdf']);
    
    %% graph properties
    % high-ranked neuron connections
    hr_nedge{n,1} = sum(pre_model.graph(hr_indx,Stim_cells),2);
    hr_nedge{n,2} = sum(post_model.graph(hr_indx,Stim_cells),2);
    
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
    
    %% plot node strength against lcc
%     nodesz = 30;
% 
%     figure; set(gcf,'color','w','position',[1983 442 339 290])
%     subplot(1,2,1); hold on
%     scatter(ep_sum{n,1},lcc{n,1},nodesz,mycc.gray_light,'filled')
%     scatter(ep_sum{n,1}(hr_indx),lcc{n,1}(hr_indx),nodesz,mycc.red,'filled')
%     xlabel('node strength'); ylabel('local clustering coeff')
%     subplot(1,2,2); hold on
%     scatter(ep_sum{n,2},lcc{n,2},nodesz,mycc.gray_light,'filled')
%     scatter(ep_sum{n,2}(hr_indx),lcc{n,2}(hr_indx),nodesz,mycc.red,'filled')
%     xlabel('node strength'); ylabel('local clustering coeff')
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_nodestrength_lcc.pdf']);
    
    %% predict with model
    % change single node activity and predict with LL
    LL_frame_pre = zeros(num_neuron,num_frame_pre,2);
    LL_frame_post = zeros(num_neuron,num_frame_post,2);
    for ii = 1:num_neuron
        for jj = 1:num_frame_pre
            frame_vec = pre_data(jj,:);
            frame_vec(ii) = 0;
            LL_frame_pre(ii,jj,1) = compute_avg_log_likelihood(pre_model.node_pot(1:num_neuron),...
                pre_model.edge_pot(1:num_neuron,1:num_neuron),pre_model.logZ,frame_vec);
            frame_vec(ii) = 1;
            LL_frame_pre(ii,jj,2) = compute_avg_log_likelihood(pre_model.node_pot(1:num_neuron),...
                pre_model.edge_pot(1:num_neuron,1:num_neuron),pre_model.logZ,frame_vec);
        end
        for jj = 1:num_frame_post
            frame_vec = post_data(jj,:);
            frame_vec(ii) = 0;
            LL_frame_post(ii,jj,1) = compute_avg_log_likelihood(post_model.node_pot(1:num_neuron),...
                post_model.edge_pot(1:num_neuron,1:num_neuron),post_model.logZ,frame_vec);
            frame_vec(ii) = 1;
            LL_frame_post(ii,jj,2) = compute_avg_log_likelihood(post_model.node_pot(1:num_neuron),...
                post_model.edge_pot(1:num_neuron,1:num_neuron),post_model.logZ,frame_vec);
        end
    end
    LL_on_pre = squeeze(LL_frame_pre(:,:,2)-LL_frame_pre(:,:,1));
    LL_on_post = squeeze(LL_frame_post(:,:,2)-LL_frame_post(:,:,1));
    
    % make predictions
    thr_pre = zeros(num_neuron,1);
    thr_post = zeros(num_neuron,1);
    LL_pred_pre = nan(num_neuron,num_frame_pre);
    LL_pred_post = nan(num_neuron,num_frame_post);
    for ii = 1:num_neuron
        [LL_pred_pre(ii,:),thr_pre(ii)] = pred_from_LL(LL_on_pre(ii,:),qnoise);
        [LL_pred_post(ii,:),thr_post(ii)] = pred_from_LL(LL_on_post(ii,:),qnoise);
    end
    
    % calculate TPR and FPR
    TPR_pre = zeros(num_neuron,1);
    FPR_pre = zeros(num_neuron,1);
    TPR_post = zeros(num_neuron,1);
    FPR_post = zeros(num_neuron,1);
    for jj = 1:num_neuron
        TP = sum(LL_pred_pre(jj,:)==1&opto_stim_pre'==1);
        FP = sum(LL_pred_pre(jj,:)~=1&opto_stim_pre'==1);
        TN = sum(LL_pred_pre(jj,:)~=1&opto_stim_pre'~=1);
        FN = sum(LL_pred_pre(jj,:)==1&opto_stim_pre'~=1);
        TPR_pre(jj) = TP/(TP+FN);
        FPR_pre(jj) = FP/(FP+TN);
        
        TP = sum(LL_pred_post(jj,:)==1&opto_stim_post'==1);
        FP = sum(LL_pred_post(jj,:)~=1&opto_stim_post'==1);
        TN = sum(LL_pred_post(jj,:)~=1&opto_stim_post'~=1);
        FN = sum(LL_pred_post(jj,:)==1&opto_stim_post'~=1);
        TPR_post(jj) = TP/(TP+FN);
        FPR_post(jj) = FP/(FP+TN);
    end
    
    % plot scatter in ROC space
    nodesz = 15;
    figure; set(gcf,'color','w','position',[2023 322 445 211])
    subplot(1,2,1);hold on
    plot([0 1],[0 1],'k--')
    scatter(FPR_pre,TPR_pre,nodesz,mycc.gray,'filled')
    scatter(FPR_pre(hr_indx),TPR_pre(hr_indx),nodesz,mycc.red,'filled')
    xlim([0 1]); ylim([0 1])
    xlabel('FPR'); ylabel('TPR'); title('pre')
    subplot(1,2,2);hold on
    plot([0 1],[0 1],'k--')
    scatter(FPR_post,TPR_post,nodesz,mycc.gray,'filled')
    scatter(FPR_post(hr_indx),TPR_post(hr_indx),nodesz,mycc.red,'filled')
    xlim([0 1]); ylim([0 1])
    xlabel('FPR'); ylabel('TPR'); title('post')
    
    saveas(gcf,[fig_path 'ROC_scatter.pdf']);
    
end

% save([save_path 'opto_spont_prop.mat'],'-v7.3');

%% distribution of P(k)
nsz = 5;
figure;
subplot(2,3,1); hold on
for ii = 1:num_expt
    plot(ndeg_bin_range,ndeg_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
    plot(ndeg_bin_range,ndeg_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
    scatter(ndeg_bin_range,ndeg_dist(ii,:,1),nsz,mycc.gray,'o','filled')
    scatter(ndeg_bin_range,ndeg_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
end

subplot(2,3,2); hold on
for ii = 1:num_expt
    plot(ndeg_bin_range,ndeg_in_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
    plot(ndeg_bin_range,ndeg_in_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
    scatter(ndeg_bin_range,ndeg_in_dist(ii,:,1),nsz,mycc.gray,'o','filled')
    scatter(ndeg_bin_range,ndeg_in_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
end

subplot(2,3,3); hold on
for ii = 1:num_expt
    plot(ndeg_bin_range,ndeg_out_dist(ii,:,1),'color',mycc.gray,'linewidth',linew);
    plot(ndeg_bin_range,ndeg_out_dist(ii,:,2),'color',mycc.blue_light,'linewidth',linew);
    scatter(ndeg_bin_range,ndeg_out_dist(ii,:,1),nsz,mycc.gray,'o','filled')
    scatter(ndeg_bin_range,ndeg_out_dist(ii,:,2),nsz,mycc.blue_light,'o','filled')
end

subplot(2,3,4); hold on
for ii = 1:num_expt
    plot(log(ndeg_bin_range),log(ndeg_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
    plot(log(ndeg_bin_range),log(ndeg_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
    scatter(log(ndeg_bin_range),log(ndeg_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
    scatter(log(ndeg_bin_range),log(ndeg_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
end

subplot(2,3,5); hold on
for ii = 1:num_expt
    plot(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
    plot(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
    scatter(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
    scatter(log(ndeg_bin_range),log(ndeg_in_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
end

subplot(2,3,6); hold on
for ii = 1:num_expt
    plot(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,1)),'color',mycc.gray,'linewidth',linew);
    plot(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,2)),'color',mycc.blue_light,'linewidth',linew);
    scatter(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,1)),nsz,mycc.gray,'o','filled')
    scatter(log(ndeg_bin_range),log(ndeg_out_dist(ii,:,2)),nsz,mycc.blue_light,'o','filled')
end

%% graph properties - whole network
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
set(gca,'xtick',[0.5 1],'xticklabel',{'pre','post'},'linewidth',linew)
set(gca,'position',gcapos);
box off

% high-ranked neuron connections
subplot(2,3,2);hold on;
pp = plot_pair_graph(cell2mat(hr_nedge(:,1)),cell2mat(hr_nedge(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('high-ranked neuron connections'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% edge pot sum
subplot(2,3,3)
pp = plot_pair_graph(cell2mat(ep_sum(:,1)),cell2mat(ep_sum(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('rank'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(2,3,4)
pp = plot_pair_graph(cell2mat(ndeg(:,1)),cell2mat(ndeg(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node degree');title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(2,3,5)
pp = plot_pair_graph(cell2mat(lcc(:,1)),cell2mat(lcc(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(2,3,6)
pp = plot_pair_graph(cell2mat(cent(:,1)),cell2mat(cent(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('all')
saveas(gcf,[fig_path 'opto_spont_graph_prop_on_unnormalized_on.pdf']);

%% stim network
figure; set(gcf,'color','w','position',[2006 320 691 425])

% density
boxwd = 0.2;
subplot(2,3,1);hold on;
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
subplot(2,3,2)
pp = plot_pair_graph(cell2mat(ep_sum_in(:,1)),cell2mat(ep_sum_in(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position'); title(num2str(pp));
ylabel('rank');
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(2,3,3)
pp = plot_pair_graph(cell2mat(ndeg_in(:,1)),cell2mat(ndeg_in(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node degree'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(2,3,4)
pp = plot_pair_graph(cell2mat(lcc_in(:,1)),cell2mat(lcc_in(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(2,3,5)
pp = plot_pair_graph(cell2mat(cent_in(:,1)),cell2mat(cent_in(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('stim')
saveas(gcf,[fig_path 'opto_spont_graph_prop_stim_network_unnormalized_on.pdf']);

%% nostim network
figure; set(gcf,'color','w','position',[2006 320 691 425])

% density
boxwd = 0.2;
subplot(2,3,1);hold on;
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
subplot(2,3,2)
pp = plot_pair_graph(cell2mat(ep_sum_out(:,1)),cell2mat(ep_sum_out(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('rank'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% node degree
subplot(2,3,3)
pp = plot_pair_graph(cell2mat(ndeg_out(:,1)),cell2mat(ndeg_out(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('node degree'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% lcc
subplot(2,3,4)
pp = plot_pair_graph(cell2mat(lcc_out(:,1)),cell2mat(lcc_out(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('clustering coeff'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

% centrality
subplot(2,3,5)
pp = plot_pair_graph(cell2mat(cent_out(:,1)),cell2mat(cent_out(:,2)),mycc.black,mycc.blue,p);
gcapos = get(gca,'position');
ylabel('centrality'); title(num2str(pp));
set(gca,'position',gcapos);
legend off; box off

suptitle('no stim')
saveas(gcf,[fig_path 'opto_spont_graph_prop_nostim_network_unnormalized_on.pdf']);


end