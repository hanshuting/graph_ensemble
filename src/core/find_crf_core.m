function [ens_crf,results] = find_crf_core(best_model,shuffle_model,data,vis_stim)
% Find CRF core neurons

%% parameters
num_stim = length(unique(setdiff(vis_stim,0)));
num_node = size(best_model.graph,1);
num_frame = length(vis_stim);
    
%% edge potential sum
% real model
best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
best_model.ep_on = best_model.ep_on + best_model.ep_on';
ns = sum(best_model.ep_on,2)/2;
ns(sum(best_model.graph,2)==0) = NaN;

% shuffled models
for ii = 1:length(shuffle_model.graphs)
    shuffle_model.ep_on{ii} = getOnEdgePot(shuffle_model.graphs{ii},...
        shuffle_model.G{ii});
    shuffle_model.ep_on{ii} = shuffle_model.ep_on{ii} + shuffle_model.ep_on{ii}';
    shuffle_model.epsum{ii} = sum(shuffle_model.ep_on{ii},2)/2;
    shuffle_model.epsum{ii}(sum(shuffle_model.graphs{ii},2)==0) = NaN;
end
ns_shuff = nanmean(cellfun(@(x) nanmean(x),shuffle_model.epsum));
ns_shuff_sd = nanstd(cellfun(@(x) nanmean(x),shuffle_model.epsum));

%% find ensemble with CRF
% predict each neuron in turn
LL_frame = zeros(num_node,num_frame,2);
for ii = 1:num_node
    for jj = 1:num_frame
        frame_vec = data(:,jj)';
        frame_vec(ii) = 0;
        LL_frame(ii,jj,1) = compute_avg_log_likelihood(best_model.node_pot,...
            best_model.edge_pot,best_model.logZ,frame_vec);
        frame_vec(ii) = 1;
        LL_frame(ii,jj,2) = compute_avg_log_likelihood(best_model.node_pot,...
            best_model.edge_pot,best_model.logZ,frame_vec);
    end
end
LL_on = squeeze(LL_frame(:,:,2)-LL_frame(:,:,1));

% calculate AUC for each node and each stim
auc = zeros(num_node,num_stim);
for ii = 1:num_stim
    true_label = double(vis_stim==ii)';
    for jj = 1:num_node
        [~,~,~,auc(jj,ii)] = perfcurve(true_label,LL_on(jj,:),1);
    end
end

% find core by AUC and node strength
auc_ctrl = cell(num_stim,1);
ens_crf = cell(num_stim,1);
for ii = 1:num_stim
    num_ens = sum(best_model.graph(num_node-num_stim+ii,:));
    for jj = 1:100
        rd_ens = randperm(num_node,num_ens);
        [~,sim_core] = core_cos_sim(rd_ens,data',...
            true_label);
        [~,~,~,auc_ctrl{ii}(jj)] = perfcurve(true_label,sim_core,1);
    end
    ens_crf{ii} = find(auc(:,ii)>(mean(auc_ctrl{ii})+std(auc_ctrl{ii}))&...
        (ns>(ns_shuff+ns_shuff_sd)));
    ens_crf{ii} = setdiff(ens_crf{ii},num_node-num_stim+ii);
end
  
%% plot AUC and node strength
nodesz = 30;
nsmi = min(ns);
nsma = max(ns);
aucmi = 0;
aucma = 1;
figure; set(gcf,'color','w','position',[1967 615 555 253])
for ii = 1:num_stim
    subplot(1,num_stim,ii); hold on
    scatter(ns,auc(:,ii),nodesz,0.5*[1 1 1],'filled')
    scatter(ns(ens_crf{ii}),auc(ens_crf{ii},ii),nodesz,[1 0.4 0.4],'filled')
%         plot([nsmi nsma],th(ii)*[1 1],'k--');
    plot([nsmi nsma],mean(auc_ctrl{ii})*[1 1],'k--');
    plot([nsmi nsma],(mean(auc_ctrl{ii})+std(auc_ctrl{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot([nsmi nsma],(mean(auc_ctrl{ii})-std(auc_ctrl{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot(ns_shuff*[1 1],[aucmi aucma],'k--');
    plot((ns_shuff+ns_shuff_sd)*[1 1],[aucmi aucma],'--',...
        'color',0.7*[1 1 1]);
    plot((ns_shuff-ns_shuff_sd)*[1 1],[aucmi aucma],'--',...
        'color',0.7*[1 1 1]);
    xlim([nsmi nsma]); ylim([aucmi aucma])
    xlabel('node strength'); ylabel(['AUC' num2str(ii)]);
end
    
%% package results
results.auc = auc;
results.auc_ctrl = auc_ctrl;
results.best_model = best_model;
results.ens_crf = ens_crf;
results.ns = ns;
results.ns_shuff = ns_shuff;
results.ns_shuff_sd = ns_shuff_sd;
results.data = logical(data);
results.LL_frame = LL_frame;
results.LL_on = LL_on;
results.shuffle_model = shuffle_model;
results.vis_stim = vis_stim;


end