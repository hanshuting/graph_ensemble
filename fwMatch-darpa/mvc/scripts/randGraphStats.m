%% --------- random graph properties --------------

% parameters
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high','on_high'}};
% ee = {{'01','02'},{'vis_01_all','vis_02_all','vis_03_all','vis_04_all'}};
data_path = 'C:\Shuting\fwMatch\data\';

num_rand = 100;
k = 4;

%% cc model

for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\'];
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\'];
        comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
        save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\rand_graph\'];
        
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        graph = cc_graph;
        
        %% generate random graphs
        num_node = size(graph,1);
        num_edge = sum(graph(:))/2;
        cc_rand_graph = zeros(num_node,num_node,num_rand);
    
        for i = 1:num_rand
            cc_rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end
    
        %% clique and community properties
        num_rand_clique = zeros(length(k),num_rand);
        num_rand_comm = zeros(length(k),num_rand);
        sz_rand_comm = cell(length(k),num_rand);
        for i = 1:length(k)
            
            fprintf('k=%u\n',k(i));
            
            clique_cc_rand = cell(num_rand,1);
            comm_cc_rand = cell(num_rand,1);
            for j = 1:num_rand
    
%                 fprintf('random graph #%u\n',j);
                [comm,clique,~] = k_clique(k(i),cc_rand_graph(:,:,j));
                clique_cc_rand{j} = clique;
                comm_cc_rand{j} = comm;
                
                % number of k cliques
                sz_clique = cellfun('length',clique);
                sz_clique = sz_clique(sz_clique>k(i));
                if ~isempty(sz_clique)
                    sz_clique = num2cell(sz_clique);
                    num_clique = cellfun(@(x) nchoosek(x,k(i)),sz_clique,'uniformoutput',false);
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
                num2str(k(i)) '.mat'],'clique_cc_rand','comm_cc_rand',...
                'cc_rand_graph','-v7.3');
            save([save_path expt_name{n} '_' expt_ee{e} '_cc_rand_graph_stats_k_' ...
                num2str(k(i)) '.mat'],'num_rand_clique','num_rand_comm','sz_rand_comm','-v7.3');
            
        end
        
        
    end

end


%% loopy model
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\'];
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\'];
        comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
        save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\rand_graph\'];
        
        load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        
        %% generate random graphs
        num_node = size(graph,1);
        num_edge = sum(graph(:))/2;
        loopy_rand_graph = zeros(num_node,num_node,num_rand);
    
        for i = 1:num_rand
            loopy_rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end
    
        %% clique and community properties
        num_rand_clique = zeros(length(k),num_rand);
        num_rand_comm = zeros(length(k),num_rand);
        sz_rand_comm = cell(length(k),num_rand);
        for i = 1:length(k)

            fprintf('k=%u\n',k(i));
            
            clique_loopy_rand = cell(num_rand,1);
            comm_loopy_rand = cell(num_rand,1);
            for j = 1:num_rand
    
%                 fprintf('random graph #%u\n',j);
                [comm,clique,~] = k_clique(k(i),loopy_rand_graph(:,:,j));
                comm_loopy_rand{j} = comm;
                clique_loopy_rand{j} = clique;
                
                % number of k cliques
                sz_clique = cellfun('length',clique);
                sz_clique = sz_clique(sz_clique>k(i));
                if ~isempty(sz_clique)
                    sz_clique = num2cell(sz_clique);
                    num_clique = cellfun(@(x) nchoosek(x,k(i)),sz_clique,'uniformoutput',false);
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
            
            save([comm_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
                num2str(k(i)) '.mat'],'clique_loopy_rand','comm_loopy_rand',...
                'loopy_rand_graph','-v7.3');
            save([save_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_stats_k_' ...
                num2str(k(i)) '.mat'],'num_rand_clique','num_rand_comm','sz_rand_comm','-v7.3');
            
        end
    
    
    end

end
