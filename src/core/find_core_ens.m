function [core_nodes, results] = find_core_ens(best_model,shuffle_model,data,stimuli, num_controls)
% FIND_TEMPORAL_ENS_NODES Find multi-timeframe neuron ensembles with CRF models.
%
% Input
%   best_model: Trained CRF model to analyze.
%   shuffle_model: Shuffled dataset CRFs.
%   data: Expected to be a timeframes by neurons binary matrix.
%   stimuli: Expected to be a timeframes by stimuli binary matrix.
%   num_controls: Optional. Number of random ensembles used to generate control
%       statistics for each stimulus.
% Output
%   ens_nodes: Cell vector containing the ensemble nodes found for each stimuli.
%       Each cell is a cell vector of length time_span containing the ensemble
%       detected at each offset frame of the time_span window.
%   results: Details used to identify the ensembles.
    if nargin < 5
        num_controls = 100;
    end

    [num_frame, num_stim] = size(stimuli);
    if all(all(data(:, end-num_stim+1:end) == stimuli))
        warning('Last %d neurons in data perfectly match stimuli. Be sure data contains only neuron recordings.', ...
                num_stim)
    end
    num_node = size(best_model.graph,1);
    num_orig_neuron = size(data, 2);
    time_span = best_model.time_span;

    % expand for additional time_span
    if time_span > 1
        orig_data = data;
        data = add_lookback_nodes(orig_data, time_span);
    end
    data = [data stimuli];
    assert(size(data,1) == num_frame, 'data frames do not match stimuli frames')
    assert(size(data,2) == num_node, 'data nodes do not match best_model.graph')

    % Ensure graph is symmetric
    if nnz(best_model.graph - logical(best_model.edge_pot)) ~= 0
        [~, best_name] = fileparts(params.best_fname);
        fprintf('Forcing graph equal to logical(edge_pot) for best_model, %s.\n', ...
            best_name);
        best_model.graph = logical(best_model.edge_pot);
    end
    
    % calculate edge potential sum
    best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
    % getOnEdgePot returns just upper triangle - copy to lower triangle to make symmetric
    best_model.ep_on = best_model.ep_on + best_model.ep_on';
    
    % assign NaN to edges that don't exist
    best_model.ep_on(best_model.graph==0) = NaN;
    
    % Sum potential for ALL edges each node is party to
    epsum = nansum(best_model.ep_on,2);
    epsum(sum(best_model.graph,2)==0) = NaN;

    % shuffled models
    for ii = 1:length(shuffle_model.graphs)
        if nnz(shuffle_model.graphs{ii} - logical(shuffle_model.edge_pot{ii})) ~= 0
            fprintf('Forcing graph equal to logical(edge_pot) for shuffle_model %d.\n', ii);
            shuffle_model.graphs{ii} = logical(shuffle_model.edge_pot{ii});
        end
        shuffle_model.ep_on{ii} = getOnEdgePot(shuffle_model.graphs{ii},...
            shuffle_model.G{ii});
        % getOnEdgePot returns just upper triangle - copy to lower triangle to make symmetric
        shuffle_model.ep_on{ii} = shuffle_model.ep_on{ii} + shuffle_model.ep_on{ii}';
        shuffle_model.epsum{ii} = nansum(shuffle_model.ep_on{ii},2);
        shuffle_model.epsum{ii}(sum(shuffle_model.graphs{ii},2)==0) = NaN;
    end
    shuffle_model.mepsum = nanmean(cellfun(@(x) nanmean(x),shuffle_model.epsum));
    shuffle_model.sdepsum = nanstd(cellfun(@(x) nanmean(x),shuffle_model.epsum));

    %% find ensembles
    ens_crf = cell(num_stim,1);
    for ii = 1:num_stim
        ens_crf{ii} = find(best_model.graph(:,num_node-num_stim+ii));
        ens_crf{ii} = setdiff(ens_crf{ii},num_node-num_stim+ii);
    end
    
    %% find core nodes in ensemble
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
        for jj = 1:num_node
            [~,~,~,auc(jj,ii)] = perfcurve(true_label(ii, :), LL_on(jj,:), 1);
        end
    end

    % calculate node strength to stimulus nodes
    ns_stim = cell(num_stim,1);
    ep_on = best_model.ep_on;
    ep_on(best_model.graph==0) = NaN;
    for ii = 1:num_stim
        ns = exp(ep_on(ens_crf{ii},ens_crf{ii}));
        ns(isnan(ns)) = 0;
        ns_stim{ii} = sum(ns,2);
    end
    
    % node strength controls
    ns_stim_ctrl = cell(num_stim,1);
    for k = 1:num_controls
        for ii = 1:num_stim
%             ens_ctrl = find(shuffle_model.graph{k}(num_node-num_stim+ii));
%             ens_ctrl = setdiff(ens_ctrl,num_node-num_stim+ii);
            ep_on = shuffle_model.ep_on{k};
            ep_on(shuffle_model.graphs{k}==0) = NaN;
            ns = exp(ep_on(ens_crf{ii},ens_crf{ii}));
            ns(isinf(ns)) = NaN;
            ns(isnan(ns)) = 0; % force 0 to be minimum given exp function
            ns_stim_ctrl{ii}(:,k) = sum(ns,2);
%             if all(isnan(ns))
%                 ns_stim_ctrl(ii,k) = NaN;
%             else
%                 ns_stim_ctrl(ii,k) = sum(ns);
%             end
        end
    end
    
    % find ensembles
    auc_rd = cell(num_stim,1);
    ns_rd = cell(num_stim,1);
    core_crf = cell(num_stim,1);
    for ii = 1:num_stim
        size_ens = length(ens_crf{ii});
        % Generate random controls for current stimulus
        ns_rd{ii} = zeros(num_controls,1);
        for jj = 1:num_controls
            rd_ens = zeros(1, num_node);
            rd_ens(randperm(length(rd_ens), size_ens)) = 1;
            % Shouldnt this only pass the population vectors from data, i.e. omit the stim nodes?
            sim_core = 1-pdist2(data,rd_ens,'cosine')';
            [~,~,~,auc_rd{ii}(jj)] = perfcurve(true_label(ii, :), sim_core, 1);
        end
        ns_thr = nanmean(ns_stim_ctrl{ii},2)+nanstd(ns_stim_ctrl{ii},[],2);
        core_crf{ii} = find(auc(ens_crf{ii},ii)>(mean(auc_rd{ii})+std(auc_rd{ii}))&...
            (ns_stim{ii})>ns_thr);
        core_crf{ii} = ens_crf{ii}(core_crf{ii});
    end

    %% Convert nodes to neurons
    core_nodes = cell(num_stim, 1);
    for ii = 1:num_stim
        core_nodes{ii} = cell(time_span, 1);
        for jj = 1:time_span
            current_nodes = (core_crf{ii} > ((jj-1) * num_orig_neuron)) & ...
                            (core_crf{ii} <= (jj * num_orig_neuron ));
            core_nodes{ii}{jj} = core_crf{ii}(current_nodes) - (jj-1) * num_orig_neuron;
        end
    end

    %% package results
    results.auc = auc;
    results.auc_ens = auc_rd;
    results.best_model = best_model;
    results.core_crf = core_crf;
    results.ens_crf = ens_crf;
    results.data = logical(data);
    results.epsum = epsum;
    results.LL_frame = LL_frame;
    results.LL_on = LL_on;
    results.shuffle_model = shuffle_model;
    results.stimuli = stimuli;
    results.time_span = time_span;

end
