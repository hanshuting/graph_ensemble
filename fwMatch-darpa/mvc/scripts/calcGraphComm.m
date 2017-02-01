function [] = calcGraphComm(param)

expt_name = param.expt_name;
ee = param.ee;
data_path = param.data_path;
result_path_base = param.result_path_base;
cc_type = param.cc_type;
comm_type = param.comm_type;
k_seq = param.k_seq;

% make xcorr graph
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % load real and random models
        model_path = [result_path_base expt_name{n} '\models\'];
        cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\'];
        save_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
        
        % load graphs
        crf_graph = load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        
        % keep all edges, or only positive edges
        crf_graph = crf_graph.graph;
        if strcmp(param.ge_meth,'on')
            num_node = size(crf_graph,1);
            num_edge = sum(sum(logical(tril(crf_graph))));
            edge_list = zeros(num_edge,2);
            [edge_list(:,2),edge_list(:,1)] = find(tril(crf_graph));
            G_crf = zeros(num_node,num_node);
            for i = 1:size(G_crf,2)
                node_1 = edge_list(i,1);
                node_2 = edge_list(i,2);
                G_crf(node_1,node_2) = G_crf(4,i);%-G(3,i)-G(2,i);
                G_crf(node_2,node_1) = G_crf(4,i);%-G(3,i)-G(2,i);
            end
            G_crf = G_crf>0;
        elseif strcmp(param.ge_meth,'full')
            G_crf = crf_graph;
        else
            error('Undefined graph edge option ''%s''\n',param.ge_meth);
        end
        
        % find communities for each k
        for k = k_seq
            
            fprintf('k=%u\n',k);
            [comm_loopy,clique_loopy,~] = k_clique(k,G_crf);
            [comm_cc,clique_cc,~] = k_clique(k,cc_graph);

            % save results
            save([save_path expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
                num2str(k) '.mat'],'comm_loopy','clique_loopy','-v7.3');
            save([save_path expt_name{n} '_' expt_ee{e} '_cc_comm_k_' ...
                num2str(k) '.mat'],'comm_cc','clique_cc','-v7.3');
            
        end
        
    end
end

end