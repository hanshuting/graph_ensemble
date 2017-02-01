function [] = randGraphStatsCC(param)

expt_name = param.expt_name;
ee = param.ee;
num_rand = param.num_rand;
data_path = param.data_path;
result_path_base = param.result_path_base;
k_seq = param.k_seq;
cc_type = param.cc_type;
comm_type = param.comm_type;
ge_meth = param.ge_meth;

%% cc model
fprintf('CC model...\n');
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\'];
        comm_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
        save_path = [result_path_base expt_name{n} '\rand_graph_' ge_meth '\'];
        
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        
        %% generate random graphs
        num_node = size(cc_graph,1);
        num_edge = sum(cc_graph(:))/2;
        cc_rand_graph = zeros(num_node,num_node,num_rand);
    
        for i = 1:num_rand
            cc_rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end
    
        %% clique and community properties
        for i = 1:length(k_seq)
            
            fprintf('k=%u\n',k_seq(i));
            
            num_rand_clique = zeros(1,num_rand);
            num_rand_comm = zeros(1,num_rand);
            sz_rand_comm = cell(1,num_rand);
            clique_cc_rand = cell(num_rand,1);
            comm_cc_rand = cell(num_rand,1);
            for j = 1:num_rand
    
%                 fprintf('random graph #%u\n',j);
                [comm,clique,~] = k_clique(k_seq(i),cc_rand_graph(:,:,j));
                clique_cc_rand{j} = clique;
                comm_cc_rand{j} = comm;
                
                % number of k cliques
                sz_clique = cellfun('length',clique);
                sz_clique = sz_clique(sz_clique>k_seq(i));
                if ~isempty(sz_clique)
                    sz_clique = num2cell(sz_clique);
                    num_clique = cellfun(@(x) nchoosek(x,k_seq(i)),sz_clique,'uniformoutput',false);
                    num_rand_clique(i,j) = sum(cell2mat(num_clique));
                else
                    num_rand_clique(i,j) = 0;
                end
    
                % number and size of communities
                if ~isempty(comm)
                    sz_rand_comm{i,j} = cellfun('length',comm);
                    num_rand_comm(i,j) = length(sz_rand_comm{i,j});
                else
                    sz_rand_comm{i,j} = 0;
                    num_rand_comm(i,j) = 0;
                end
    
            end
            
            save([comm_path expt_name{n} '_' expt_ee{e} '_cc_rand_graph_comm_k_' ...
                num2str(k_seq(i)) '.mat'],'clique_cc_rand','comm_cc_rand',...
                'cc_rand_graph','-v7.3');
            save([save_path expt_name{n} '_' expt_ee{e} '_cc_rand_graph_stats_k_' ...
                num2str(k_seq(i)) '.mat'],'num_rand_clique','num_rand_comm','sz_rand_comm','-v7.3');
            
        end 
    end
end

end
