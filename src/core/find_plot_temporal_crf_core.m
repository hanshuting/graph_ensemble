function [results] = find_plot_temporal_crf_core(best_model,shuffle_model,data,stimuli, coords)
% FIND_PLOT_TEMPORAL_CRF_CORE Find and plot multi-timeframe neuron ensembles with CRF models.
%
% Input
%   best_model: Trained CRF model to analyze.
%   shuffle_model: Shuffled dataset CRFs.
%   data: Expected to be a timeframes by neurons binary matrix.
%   stimuli: Expected to be a timeframes by stimuli binary matrix.
%   coords: Optional. 2D spatial coordinates to display each neuron.
%       Neurons by coordinates.
    if nargin < 5
        coords = [];
        num_subplots = 1;
    else
        num_subplots = 2;
    end
    % Set number of random ensembles used to generate control statistics for
    % each stimulus
    num_controls = 100;

    [num_frame, num_stim] = size(stimuli);
    num_node = size(best_model.graph,1);
    num_orig_neuron = size(data, 2);
    time_span = best_model.time_span;

    % expand for additional time_span nodes
    if time_span > 1
        coords = repmat(coords, time_span, 1);
    end

    [~, results] = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli, num_controls);
    core_crf = results.core_crf;
    epsum = results.epsum;
    auc = results.auc;
    auc_ens = results.auc_ens;
    shuffle_model = results.shuffle_model;

    %% plot
    nodesz = 30;
    nsmi = min(epsum);
    nsma = max(epsum);
    aucmi = 0;
    aucma = 1;
    f = figure; set(gcf,'color','w')
    f.Name = sprintf('K=%d', time_span);
    color_by_offset = @(x) floor((x-1)/num_orig_neuron) / max(1, time_span-1);
    for ii = 1:num_stim

        % AUC - node strength plot
        cur_axes = subplot(num_subplots,num_stim,ii); hold on
        colormap(cur_axes, autumn)
        scatter(epsum,auc(:,ii),nodesz,0.5*[1 1 1],'filled')
        % Stimuli nodes blue
        scatter(epsum(end - num_stim + 1:end),auc(end - num_stim + 1:end,ii),nodesz,[0 0 1],'filled')
        % Core nodes colored red->yellow according to how frame-offset
        scatter(epsum(core_crf{ii}),auc(core_crf{ii},ii),nodesz,arrayfun(color_by_offset, core_crf{ii}),'filled')
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
        if num_subplots > 1
            subplot(num_subplots,num_stim,ii+num_stim);
            plotGraphHighlight(coords,mod(core_crf{ii}-1, num_orig_neuron)+1,'red',1 / time_span)
        end

    end
end
