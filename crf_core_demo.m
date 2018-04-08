% DEMO SCRIPT FOR FINDING ENSEMBLE CORE NEURONS WITH CRF MODELS

%% load data, pre-process
% load data
load('data_demo.mat');
best_model = load('model_demo.mat');
shuffle_model = load('model_shuffled_demo.mat');

% extract some parameters
num_stim = length(unique(setdiff(vis_stim,0)));
num_node = size(best_model.graph,1);
num_frame = length(vis_stim);
data = data';

% calculate edge potential sum
best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
best_model.ep_on = best_model.ep_on - tril(best_model.ep_on);
epsum = sum(best_model.ep_on,2);
epsum(sum(best_model.graph,2)==0) = NaN;

% shuffled models
for ii = 1:length(shuffle_model.graphs)
    shuffle_model.ep_on{ii} = getOnEdgePot(shuffle_model.graphs{ii},...
        shuffle_model.G{ii});
    shuffle_model.epsum{ii} = sum(shuffle_model.ep_on{ii},2);
    shuffle_model.epsum{ii}(sum(shuffle_model.graphs{ii},2)==0) = NaN;
end
shuffle_model.mepsum = nanmean(cellfun(@(x) nanmean(x),shuffle_model.epsum));
shuffle_model.sdepsum = nanstd(cellfun(@(x) nanmean(x),shuffle_model.epsum));

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

% calculate AUC
true_label = zeros(num_stim, num_frame);
auc = zeros(num_node,num_stim);
for ii = 1:num_stim
    true_label(ii, :) = double(vis_stim==ii)';
    for jj = 1:num_node
        [~,~,~,auc(jj,ii)] = perfcurve(true_label(ii, :), LL_on(jj,:), 1);
    end
end

% find ensembles
auc_ens = cell(num_stim,1);
core_crf = cell(num_stim,1);
for ii = 1:num_stim
    num_ens = sum(best_model.graph(num_node-num_stim+ii,:));
    for jj = 1:100
        rd_ens = randperm(num_node,num_ens);
        [~,sim_core] = core_cos_sim(rd_ens,data', true_label(ii, :));
        [~,~,~,auc_ens{ii}(jj)] = perfcurve(true_label(ii, :), sim_core, 1);
    end
    core_crf{ii} = find(auc(:,ii)>(mean(auc_ens{ii})+std(auc_ens{ii}))&...
        (epsum>(shuffle_model.mepsum+shuffle_model.sdepsum)));
    core_crf{ii} = setdiff(core_crf{ii},num_node-num_stim+ii);
end

%% plot
nodesz = 30;
nsmi = min(epsum);
nsma = max(epsum);
aucmi = 0;
aucma = 1;
figure; set(gcf,'color','w')
for ii = 1:num_stim
    
    % AUC - node strength plot
    subplot(2,num_stim,ii); hold on
    scatter(epsum,auc(:,ii),nodesz,0.5*[1 1 1],'filled')
    % Core nodes red
    scatter(epsum(core_crf{ii}),auc(core_crf{ii},ii),nodesz,[1 0 0],'filled')
    % Stimuli nodes blue
    scatter(epsum(end - num_stim + 1:end),auc(end - num_stim + 1:end,ii),nodesz,[0 0 1],'filled')
    % Active stimulus node green
    scatter(epsum(num_node - num_stim + ii),auc(num_node - num_stim + ii,ii),nodesz,[0 1 0],'filled')
    plot([nsmi nsma],mean(auc_ens{ii})*[1 1],'k--');
    plot([nsmi nsma],(mean(auc_ens{ii})+std(auc_ens{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot([nsmi nsma],(mean(auc_ens{ii})-std(auc_ens{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot(shuffle_model.mepsum*[1 1],[aucmi aucma],'k--');
    plot((shuffle_model.mepsum+shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',0.7*[1 1 1]);
    plot((shuffle_model.mepsum-shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
        'color',0.7*[1 1 1]);
    xlim([nsmi nsma]); ylim([aucmi aucma])
    xlabel('node strength'); ylabel(['AUC' num2str(ii)]);
    title(['core #' num2str(ii)])
    
    % plot coordinates
    subplot(2,num_stim,ii+num_stim);
    plotGraphHighlight(Coord_active,core_crf{ii},[1 0 0])
    
end



