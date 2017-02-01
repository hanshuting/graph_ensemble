function [] = randGraphStatsCRF(param)

expt_name = param.expt_name;
ee = param.ee;
num_rand = param.num_rand;
data_path = param.data_path;
result_path_base = param.result_path_base;
k_seq = param.k_seq;
cc_type = param.cc_type;
comm_type = param.comm_type;
ge_meth = param.ge_meth;

% loopy model
fprintf('loopy model...\n');
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        model_path = [result_path_base expt_name{n} '\models\'];
        comm_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
        rand_graph_path = [result_path_base expt_name{n} '\rand_graph_' ge_meth '\'];
        
        crf_graph = load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        
        %% generate random graphs
        % keep all edges, or only positive edges
        crf_graph = crf_graph.graph;
        if strcmp(param.ge_meth,'on')
            num_node = size(crf_graph,1);
            num_edge = sum(sum(logical(tril(crf_graph))));
            edge_list = zeros(num_edge,2);
            [edge_list(:,2),edge_list(:,1)] = find(tril(crf_graph));
            G_crf = zeros(num_node,num_node);
            for ii = 1:size(G_crf,2)
                node_1 = edge_list(ii,1);
                node_2 = edge_list(ii,2);
                G_crf(node_1,node_2) = G_crf(4,ii);%-G(3,i)-G(2,i);
                G_crf(node_2,node_1) = G_crf(4,ii);%-G(3,i)-G(2,i);
            end
            G_crf = G_crf>0;
        elseif strcmp(param.ge_meth,'full')
            G_crf = crf_graph;
        else
            error('Undefined graph edge option ''%s''\n',param.ge_meth);
        end
        
        num_node = size(G_crf,1);
        num_edge = sum(G_crf(:))/2;
        loopy_rand_graph = zeros(num_node,num_node,num_rand);
    
        for ii = 1:num_rand
            loopy_rand_graph(:,:,ii) = mkRandomGraph(num_node,num_edge);
        end
    
        %% clique and community properties
        for ii = 1:length(k_seq)

            fprintf('k=%u\n',k_seq(ii));
            
            num_rand_clique = zeros(1,num_rand);
            num_rand_comm = zeros(1,num_rand);
            sz_rand_comm = cell(1,num_rand);
            clique_loopy_rand = cell(num_rand,1);
            comm_loopy_rand = cell(num_rand,1);
            for jj = 1:num_rand
    
%                 fprintf('random graph #%u\n',j);
                [comm,clique,~] = k_clique(k_seq(ii),loopy_rand_graph(:,:,jj));
                comm_loopy_rand{jj} = comm;
                clique_loopy_rand{jj} = clique;
                
                % number of k cliques
                sz_clique = cellfun('length',clique);
                sz_clique = sz_clique(sz_clique>k_seq(ii));
                if ~isempty(sz_clique)
                    sz_clique = num2cell(sz_clique);
                    num_clique = cellfun(@(x) nchoosek(x,k_seq(ii)),sz_clique,'uniformoutput',false);
                    num_rand_clique(ii,jj) = sum(cell2mat(num_clique));
                else
                    num_rand_clique(ii,jj) = 0;
                end
    
                % number and size of communities
                if ~isempty(comm)
                    sz_rand_comm{ii,jj} = cellfun('length',comm);
                    num_rand_comm(ii,jj) = length(sz_rand_comm{ii,jj});
                else
                    sz_rand_comm{ii,jj} = 0;
                    num_rand_comm(ii,jj) = 0;
                end
    
            end
            
            save([comm_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
                num2str(k_seq(ii)) '.mat'],'clique_loopy_rand','comm_loopy_rand',...
                'loopy_rand_graph','-v7.3');
            save([rand_graph_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_stats_k_' ...
                num2str(k_seq(ii)) '.mat'],'num_rand_clique','num_rand_comm','sz_rand_comm','-v7.3');
            
        end
    
    
    end

end

