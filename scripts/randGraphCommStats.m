function [] = randGraphCommStats(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
k_seq = param.k_seq;
% cc_type = param.cc_type;
% comm_type = param.comm_type;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.cc_path;

load(ccode_path);

%% load data
expt_count = 0;
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        expt_count = expt_count+1;
        fprintf('Processing %s_%s...\n',expt_name{n},expt_ee{e});
        
         % set path
        model_path = [result_path_base expt_name{n} '\models\']; 
%         cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\']; 
%         comm_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
        
        % load data
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
%         load([cc_path  expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
%         load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
%             num2str(k) '.mat']); % clique_loopy, comm_loopy
%         load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_comm_k_' ...
%             num2str(k) '.mat']); % clique_cc, comm_cc
%         load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_rand_graph_comm_k_' ...
%             num2str(k) '.mat']); % clique_cc_rand, comm_cc_rand
%         load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
%             num2str(k) '.mat']); % clique_loopy_rand, comm_loopy_rand
        
        %% generate random graphs
        num_node = size(graph,1);
        num_edge = sum(graph(:))/2;
        rand_graph = zeros(num_node,num_node,num_rand);

        for i = 1:num_rand
            rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end

        %% clique and community properties
        num_rand_clique = zeros(length(k),num_rand);
        num_rand_comm = zeros(length(k),num_rand);
        sz_rand_comm = cell(length(k),num_rand);
        for i = 1:length(k_seq)
            fprintf('k=%u\n;',k_seq(i));
    %         % number of node that can form a complete graph
    %         Nc = floor(sqrt(2*num_edge));
    % 
    %         % number of extra node that can form k cliques
    %         Ne = floor((num_edge-Nc)/k(i));
    %         
    %         % maximum number of k cliques
    %         Nk = nchoosek(Nc,k(i))+nchoosek(Ne*k(i),k(i)-1);

            for j = 1:num_rand

                fprintf('random graph #%u\n',j);
                [comm,clique,~] = k_clique(k_seq(i),rand_graph(:,:,j));

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
        end

        save([save_path ee{n} '_rand_graph_stats.mat'],'num_rand_clique',...
            'num_rand_comm','sz_rand_comm','-v7.3');

    end
    
end
