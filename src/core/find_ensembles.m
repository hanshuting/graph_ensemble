function [results] = find_ensembles(best_model,shuffle_model,data,stimuli, num_controls)
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

    % Sum potential for ALL edges each node is party to - exp and sum
    best_model.ep_on(best_model.graph==0) = NaN;
    ns = exp(best_model.ep_on);
    ns(isinf(ns)) = NaN;
    ns(isnan(ns)) = 0; % force 0 to be minimum given exp function
    ns = sum(ns,2);

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
        shuffle_model.ns{ii} = nansum(shuffle_model.ep_on{ii},2);
        shuffle_model.ns{ii}(sum(shuffle_model.graphs{ii},2)==0) = NaN;
        shuffle_model.ep_on{ii}(shuffle_model.graphs{ii}==0) = NaN;
        shuffle_model.ns{ii} = exp(shuffle_model.ep_on{ii});
        shuffle_model.ns{ii}(isinf(shuffle_model.ns{ii})) = NaN;
        shuffle_model.ns{ii}(isnan(shuffle_model.ns{ii})) = 0; % force 0 to be minimum given exp function
        shuffle_model.ns{ii} = sum(shuffle_model.ns{ii},2);
    end
    shuffle_model.mns = nanmean(cellfun(@(x) nanmean(x),shuffle_model.ns));
    shuffle_model.sdns = nanstd(cellfun(@(x) nanmean(x),shuffle_model.ns));

    % node potential controls
    ns_shuff = cell2mat(shuffle_model.ns);
        
    %% AUC criteria
    % predict each neuron in turn
    LL_frame = zeros(num_node,num_frame,2);
    edge_pot = best_model.edge_pot - tril(best_model.edge_pot);
    vals = [0,1];
    
    for ii = 1:num_node
        
        data_mat = data;
        data_mat(ii,:) = 1;
        LL_frame(ii,:,2) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,data_mat);
        
        for jj = 1:num_frame
            
            frame_vec = data(jj,:);
            
            for kk = 1:2
                frame_vec(ii) = vals(kk);

                % calculate node and edge effects
                log_likelihood = sum(frame_vec*best_model.node_pot)+...
                    sum(sum(edge_pot.*(frame_vec'*frame_vec)));

                % take the average
                LL_frame(ii,jj,kk) = log_likelihood-best_model.logZ;
            end
            
        end
    end
%     LL_on = squeeze(LL_frame(:,:,2)-LL_frame(:,:,1));
    LL_on = exp(squeeze(LL_frame(:,:,2)-LL_frame(:,:,1)));
    
    % calculate AUC
    true_label = stimuli';
    auc = zeros(num_node,num_stim);
    for ii = 1:num_stim
        for jj = 1:num_node
            [~,~,~,auc(jj,ii)] = perfcurve(true_label(ii, :), LL_on(jj,:), 1);
        end
    end
    
    % AUC controls
    num_ctrl_subset = 10; % use a subset due to the long runtime
    subset_ind = randperm(num_controls,num_ctrl_subset);
    auc_ctrl = zeros(num_node,num_stim,num_ctrl_subset);
    for k = 1:num_ctrl_subset
        
        indx = subset_ind(k);
        edge_pot_shuff = shuffle_model.edge_pot{indx} - tril(shuffle_model.edge_pot{indx});
        
        for ii = 1:num_node
            
            % calculate frame likelihood
            LL_frame_ctrl = zeros(num_frame,2);
            for jj = 1:num_frame
                frame_vec = data(jj,:);
                for kk = 1:2
                    frame_vec(ii) = vals(kk);
                    log_likelihood = sum(frame_vec*shuffle_model.node_pot{indx})+...
                        sum(sum(edge_pot_shuff.*(frame_vec'*frame_vec)));
                    LL_frame_ctrl(jj,kk) = log_likelihood-shuffle_model.logZ(indx);
                end
            end
%             LL_on_ctrl = squeeze(LL_frame_ctrl(:,2)-LL_frame_ctrl(:,1));
            LL_on_ctrl = exp(squeeze(LL_frame_ctrl(:,2)-LL_frame_ctrl(:,1)));
            
            % compute AUC
            for ss = 1:num_stim
                [~,~,~,auc_ctrl(ii,ss,k)] = perfcurve(true_label(ss,:),LL_on_ctrl,1);
            end
            
        end
    end
    
    %% define ensembles
    ens_crf = cell(num_stim,1);
    for ii = 1:num_stim
        
%         ens_crf{ii} = find((ns>quantile(ns_shuff,0.95,2)) & ...
%             (auc(:,ii)>quantile(squeeze(auc_ctrl(:,ii,:)),0.95,2)));
%         ens_crf{ii} = find(best_model.graph(end-num_stim+ii,:));

%         ens_crf{ii} = find((auc(:,ii)>quantile(squeeze(auc_ctrl(:,ii,:)),0.7,2))&...
%             (ns>quantile(ns_shuff,0.7,2)));
        
%         ens_crf{ii} = find((auc(:,ii)>(mean(squeeze(auc_ctrl(:,ii,:)),2)+...
%             std(squeeze(auc_ctrl(:,ii,:)),[],2)))&...
%             (ns>(mean(ns_shuff,2)+std(ns_shuff,[],2))));
        
        ens_crf{ii} = find(auc(:,ii)>mean(squeeze(auc_ctrl(:,ii,:)),2) & ...
            ns>mean(ns_shuff,2));
        
        ens_crf{ii} = setdiff(ens_crf{ii},(num_node-num_stim+1:num_node));
        
    end
    
    %% find core nodes in ensemble
    % calculate node strength within ensembles
    ns_ens = cell(num_stim,1);
    for ii = 1:num_stim
        % add up all edge potential effect for coactivity level
         ns_k = nansum(best_model.ep_on(ens_crf{ii},ens_crf{ii}),2);
         ns_k(ns_k==0) = NaN;
         ns_k = exp(ns_k); ns_k(isnan(ns_k)) = 0;
         % add node potential effect to reflect spontaneous activity level
         np = best_model.node_pot(ens_crf{ii}); np(np==0) = NaN; np = exp(np); np(isnan(np)) = 0;
         ns_ens{ii} = ns_k.*np;
%         ns_k = exp(best_model.ep_on(ens_crf{ii},ens_crf{ii}));
%         ns_ens{ii}(isinf(ns_ens{ii})) = NaN;
%         ns_ens{ii}(isnan(ns_ens{ii})) = 0;
%         ns_ens{ii} = sum(ns_ens{ii},2);
    end
    
    % node strength controls with stimulus nodes, and within ensemble
    ns_ctrl = cell(num_stim,1);
    for ii = 1:num_stim
        for k = 1:num_controls
             ns_k = nansum(shuffle_model.ep_on{k}(ens_crf{ii},ens_crf{ii}),2);
             ns_k(ns_k==0) = NaN;
             ns_k = exp(ns_k); ns_k(isnan(ns_k)) = 0;
             np = shuffle_model.node_pot{k}(ens_crf{ii}); np(np==0) = NaN; np = exp(np); np(isnan(np)) = 0;
             ns_ctrl{ii}(:,k) = ns_k.*np;
%             ns_k = exp(shuffle_model.ep_on{k}(ens_crf{ii},ens_crf{ii}));
%             ns_k(isinf(ns_k)) = NaN;
%             ns_k(isnan(ns_k)) = 0;
%             ns_ctrl{ii}(:,k) = sum(ns_k,2);
        end
    end
    
    % find core
    core_crf = cell(num_stim,1);
    for ii = 1:num_stim
        idx = (ns_ens{ii}>mean(ns_ens{ii})) & ...
            (auc(ens_crf{ii},ii)> mean(auc(ens_crf{ii},ii)));
%         idx = (ns_ens{ii}>(mean(ns_ctrl{ii},2)+std(ns_ctrl{ii},[],2))) & ...
%             (auc(ens_crf{ii},ii)>quantile(squeeze(auc_ctrl(ens_crf{ii},ii,:)),0.9,2));
%         idx = (ns_ens{ii}>(mean(ns_ctrl{ii},2))) & ...
%             (auc(ens_crf{ii},ii)>quantile(squeeze(auc_ctrl(ens_crf{ii},ii,:)),0.9,2));
%         idx = (ns_ens{ii}>(mean(ns_ctrl{ii},2))+std(ns_ctrl{ii},[],2)) & ...
%             (auc(ens_crf{ii},ii)>(mean(squeeze(auc_ctrl(ens_crf{ii},ii,:)),2)+...
%             std(squeeze(auc_ctrl(ens_crf{ii},ii,:)),[],2)));
        core_crf{ii} = ens_crf{ii}(idx);
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
    results.auc_ctrl = auc_ctrl;
    results.best_model = best_model;
    results.core_crf = core_crf;
    results.ens_crf = ens_crf;
    results.ns = ns;
    results.ns_shuff = ns_shuff;
    results.ns_ens = ns_ens;
    results.ns_ens_ctrl = ns_ctrl;
    results.data = logical(data);
    results.LL_frame = LL_frame;
    results.LL_on = LL_on;
    results.LL_on_ctrl = LL_on_ctrl;
    results.shuffle_model = shuffle_model;
    results.stimuli = stimuli;
    results.time_span = time_span;

end
