function [] = threshCRFgraphs(param)
% threshold CRF graphs by edge potential

expt_name = param.expt_name;
ee = param.ee;
data_path = param.data_path;
result_path_base = param.result_path_base;

for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set path
        model_path = [result_path_base '\' expt_name{n} '\models\'];
        save_path = [result_path_base '\' expt_name{n} '\models\'];
        
        % load graphs
        crf_graph = load([model_path expt_name{n} '_' expt_ee{e} ...
            '_loopy_best_model_full.mat']);
        F = crf_graph.F;
        G = crf_graph.G;
        logZ = crf_graph.logZ;
        node_pot = crf_graph.node_pot;
        edge_pot = crf_graph.edge_pot;
        crf_graph = crf_graph.graph;
        
        % load graphs from shuffled data
        th = load([model_path 'shuffled_' expt_name{n} '_' expt_ee{e} '_loopy.mat']);
        
        % threshold
        indx = ~isnan(th.edge_pot_shuffled);
        ccdist = fitdist(exp(th.edge_pot_shuffled(indx)),'normal');
        edge_thresh = log(icdf(ccdist,param.epth));
        edge_pot = edge_pot.*(edge_pot>edge_thresh);
        graph = crf_graph==1&edge_pot>edge_thresh;
        
        num_elm = full(sum(sum(crf_graph==1&edge_pot<=edge_thresh)));
        num_ttl = full(sum(crf_graph(:)));
        fprintf('threshold is %4.2f; eliminated %4.2f%%(%u/%u) edges\n',...
            edge_thresh,num_elm/num_ttl,num_elm,num_ttl);
        
%         % change G correspondingly
%         edge_list = zeros(sum(crf_graph(:))/2,2);
%         [edge_list(:,2),edge_list(:,1)] = find(tril(crf_graph));
%         edge_list_th = zeros(sum(graph(:))/2,2);
%         [edge_list_th(:,2),edge_list_th(:,1)] = find(tril(graph));
%         for ii = 1:size(edge_list,1)
%             if ~ismember(edge_list(ii,:),edge_list_th,'rows')
%                 G(:,ii) = 0;
%             end
%         end
        
        % save results
        save([save_path expt_name{n} '_' expt_ee{e} '_loopy_best_model_thresh.mat'],...
            'graph','F','G','edge_pot','node_pot','logZ','-v7.3');
        
    end
end

end