function [results] = temporal_crf_analysis(params)
% TEMPORAL_CRF_ANALYSIS Find multi-timeframe neuron ensembles with CRF models.
%
% Input
%   Params: A struct that contains the following fields
%       data_fname: Filepath to the data .mat file; expect variable named
%           'data', `stimuli`, `Coord_active`.
%       best_fname: Filepath to the trained CRF model to analyze; expect
%           variable named 'best_model'.
%       shuffle_fname: Filepath to the shuffled dataset CRFs; expect
%           variable named 'shuffle_model'.
%       time_span: Optional. Used if no best_model.time_span does not
%           exist.
%       num_controls: Optional. Number of random ensembles used to generate
%           control statistics for each stimulus.
if nargin < 1
    params.time_span = 2;
    params.data_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
        "highActivityAnalysis.mat";
    params.best_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
        "alejandro_highActivityData-thresh-" + int2str(params.time_span) + "_loopy/results" +  ...
        "/alejandro_highActivityData-thresh-" + int2str(params.time_span) + "_loopy_best_model_full.mat";
%         "alejandro_highActivityData-stim_loopy/results-k" + int2str(params.time_span) + ...
    params.shuffle_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
        "shuffled_alejandro_highActivityData-thresh-" + int2str(params.time_span) + "_loopy/results/" + ...
        "/shuffled_alejandro_highActivityData-thresh-" + int2str(params.time_span) + "_loopy_fulldata.mat";
%     params.data_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
%         "fullAnalysis.mat";
%     params.best_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
%         "alejandro_fullData-thresh-" + int2str(params.time_span) + "_loopy/results" +  ...
%         "/alejandro_fullData-thresh-" + int2str(params.time_span) + "_loopy_best_model_full.mat";
% %         "alejandro_highActivityData-stim_loopy/results-k" + int2str(params.time_span) + ...
%     params.shuffle_fname = "/Users/jonathanshor/GitHub/graph_ensemble/data/alejandro/" + ...
%         "shuffled_alejandro_fullData-thresh-" + int2str(params.time_span) + "_loopy/results/" + ...
%         "/shuffled_alejandro_fullData-thresh-" + int2str(params.time_span) + "_loopy_fulldata.mat";
end
    %% load data, pre-process
    % load data
    load(params.data_fname);
    best_model = load(params.best_fname);
    shuffle_model = load(params.shuffle_fname);

    % extract some parameters
    if exist('params.num_controls', 'var') == 1
        num_controls = params.num_controls;
    else
        num_controls = 100;
    end
    if exist('stimuli', 'var') ~= 1
        stimuli = zeros(length(vis_stim), length(unique(setdiff(vis_stim,0))));
        for ii = 1:size(stimuli, 2)
            stimuli(:, ii) = double(vis_stim == ii);
        end
        assert(all(all(data(:, end-size(stimuli,2)+1:end) == stimuli)))
        data = data(:, 1:end-size(stimuli,2));
    end
    num_stim = size(stimuli, 2);
    num_node = size(best_model.graph,1);
    num_orig_node = size(Coord_active, 1);
    num_frame = size(stimuli, 1);
    if exist('best_model.time_span', 'var') == 1
        time_span = best_model.time_span;
    else
        time_span = params.time_span;
    end

    % expand for additional time_span
    if time_span > 1
        orig_data = data;
        data = add_lookback_nodes(orig_data, orig_data, time_span);
        Coord_active = repmat(Coord_active, time_span, 1);
    end
    data = [data stimuli];
    assert(size(data,1) == num_frame, "data frames do not match stimuli frames")
    assert(size(data,2) == num_node, "data nodes do not match best_model.graph")

    % calculate edge potential sum
    if nnz(best_model.graph - logical(best_model.edge_pot)) ~= 0
        [~, best_name] = fileparts(params.best_fname);
        fprintf("Forcing graph equal to logical(edge_pot) for best_model, %s.\n", ...
            best_name);
        best_model.graph = logical(best_model.edge_pot);
    end
    best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
    % getOnEdgePot returns just upper triangle - need to restore symmetry
    best_model.ep_on = best_model.ep_on + best_model.ep_on';
    % In order to sum potential  for ALL edges each node is party to
    epsum = sum(best_model.ep_on,2);
    epsum(sum(best_model.graph,2)==0) = NaN;

    % shuffled models
    for ii = 1:length(shuffle_model.graphs)
        if nnz(shuffle_model.graphs{ii} - logical(shuffle_model.edge_pot{ii})) ~= 0
            [~, shuffle_name] = fileparts(params.shuffle_fname);
            fprintf("Forcing graph equal to logical(edge_pot) for shuffle_model %d, %s.\n", ...
                ii, shuffle_name);
            shuffle_model.graphs{ii} = logical(shuffle_model.edge_pot{ii});
        end
        shuffle_model.ep_on{ii} = getOnEdgePot(shuffle_model.graphs{ii},...
            shuffle_model.G{ii});
        shuffle_model.ep_on{ii} = shuffle_model.ep_on{ii} + shuffle_model.ep_on{ii}';
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
            frame_vec = data(jj,:);
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
    true_label = stimuli';
    auc = zeros(num_node,num_stim);
    for ii = 1:num_stim
%         true_label(ii, :) = double(vis_stim==ii)';
        for jj = 1:num_node
            [~,~,~,auc(jj,ii)] = perfcurve(true_label(ii, :), LL_on(jj,:), 1);
        end
    end

    % find ensembles
    auc_ens = cell(num_stim,1);
    core_crf = cell(num_stim,1);
    for ii = 1:num_stim
        size_ens = sum(best_model.graph(num_node-num_stim+ii,:));
        % Generate random controls for current stimulus
        for jj = 1:num_controls
            rd_ens = zeros(1, num_node);
            rd_ens(randperm(length(rd_ens), size_ens)) = 1;
%             [~,sim_core] = core_cos_sim(rd_ens,data, true_label(ii, :));
            % Shouldnt this only pass the population vectors from data, i.e. omit the stim nodes?
            sim_core = 1-pdist2(data,rd_ens,'cosine')';
            [~,~,~,auc_ens{ii}(jj)] = perfcurve(true_label(ii, :), sim_core, 1);
        end
        core_crf{ii} = find(auc(:,ii)>(mean(auc_ens{ii})+std(auc_ens{ii}))&...
            (epsum>(shuffle_model.mepsum+shuffle_model.sdepsum)));
        core_crf{ii} = setdiff(core_crf{ii},num_node-num_stim+1:num_node);
    end

    %% package results
    results.auc = auc;
    results.auc_ens = auc_ens;
    results.best_model = best_model;
    results.core_crf = core_crf;
    results.data = logical(data);
    results.epsum = epsum;
    results.LL_frame = LL_frame;
    results.LL_on = LL_on;
    results.shuffle_model = shuffle_model;
    results.stimuli = stimuli;
    results.time_span = time_span;

    %% plot
    nodesz = 30;
    nsmi = min(epsum);
    nsma = max(epsum);
    aucmi = 0;
    aucma = 1;
    f = figure; set(gcf,'color','w')
    [~, data_title] = fileparts(params.data_fname);
    f.Name = data_title + " K=" + time_span;
    color_by_offset = @(x) floor((x-1)/num_orig_node) / max(1, time_span-1);
    for ii = 1:num_stim

        % AUC - node strength plot
        cur_axes = subplot(2,num_stim,ii); hold on
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
        subplot(2,num_stim,ii+num_stim);
        plotGraphHighlight(Coord_active,mod(core_crf{ii}-1, num_orig_node)+1,'red',1 / time_span)

    end
end
