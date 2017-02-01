
% parameters
rng(1000);
expt_name = {'m21_d2_vis','mf2_v1'};
ee = {{'01','02'},{'vis_01_all','vis_02_all','vis_04_all'}};
num_rand = 100;
data_path = 'C:\Shuting\fwMatch\data\';
k = 3;
p = 0.05;

% calculate community size threshold
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set save path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\'];
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\'];
        save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
        
        %---------- loopy data threshold ----------%
        load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        num_node = size(graph,1);
        num_edge = sum(graph(:))/2;
        rand_graph = zeros(num_node,num_node,num_rand);
        
        % make random graph
        for i = 1:num_rand
            rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end
        
        % calculate communities
         num_rand_clique = zeros(num_rand,1);
        num_rand_comm = zeros(num_rand,1);
        sz_rand_comm = cell(num_rand,1);
    
        for i = 1:num_rand
        
            fprintf('random graph #%u\n',i);
            [comm,clique,~] = k_clique(k,rand_graph(:,:,i));
        
            % number of k cliques
            sz_clique = cellfun('length',clique);
            sz_clique = sz_clique(sz_clique>k);
            if ~isempty(sz_clique)
                sz_clique = num2cell(sz_clique);
                num_clique = cellfun(@(x) nchoosek(x,k),sz_clique,'uniformoutput',false);
                num_rand_clique(i) = sum(cell2mat(num_clique));
            else
                num_rand_clique(i) = 0;
            end
        
            % number and size of communities
            if ~isempty(comm)
                sz_rand_comm{i} = cellfun('length',comm);
                num_rand_comm(i) = length(sz_rand_comm{i});
            else
                sz_rand_comm{i} = 0;
                num_rand_comm(i) = 0;
            end
    
        end
        
        % calculate threshold
        comm_sz_mat = cell2mat(sz_rand_comm);
        bin_range = 0:1:max(comm_sz_mat);
        hist_comm = hist(comm_sz_mat,bin_range);
        hist_comm = hist_comm/sum(hist_comm);
        cdf_comm_rand = cumsum(hist_comm);
        comm_sz_thresh = find(cdf_comm_rand>1-p,1);
        
        % save results
        save([save_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
            num2str(k) '.mat'],'num_rand_clique','num_rand_comm',...
            'sz_rand_comm','comm_sz_thresh','-v7.3');
        
        %---------- cc data threshold ----------%
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        graph = cc_graph;
        num_node = size(graph,1);
        num_edge = sum(graph(:))/2;
        rand_graph = zeros(num_node,num_node,num_rand);
        
        % make random graph
        for i = 1:num_rand
            rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
        end
        
        % calculate communities
         num_rand_clique = zeros(num_rand,1);
        num_rand_comm = zeros(num_rand,1);
        sz_rand_comm = cell(num_rand,1);
    
        for i = 1:num_rand
        
            fprintf('random graph #%u\n',i);
            [comm,clique,~] = k_clique(k,rand_graph(:,:,i));
        
            % number of k cliques
            sz_clique = cellfun('length',clique);
            sz_clique = sz_clique(sz_clique>k);
            if ~isempty(sz_clique)
                sz_clique = num2cell(sz_clique);
                num_clique = cellfun(@(x) nchoosek(x,k),sz_clique,'uniformoutput',false);
                num_rand_clique(i) = sum(cell2mat(num_clique));
            else
                num_rand_clique(i) = 0;
            end
        
            % number and size of communities
            if ~isempty(comm)
                sz_rand_comm{i} = cellfun('length',comm);
                num_rand_comm(i) = length(sz_rand_comm{i});
            else
                sz_rand_comm{i} = 0;
                num_rand_comm(i) = 0;
            end
    
        end
        
        % calculate threshold
        comm_sz_mat = cell2mat(sz_rand_comm);
        bin_range = 0:1:max(comm_sz_mat);
        hist_comm = hist(comm_sz_mat,bin_range);
        hist_comm = hist_comm/sum(hist_comm);
        cdf_comm_rand = cumsum(hist_comm);
        comm_sz_thresh = find(cdf_comm_rand>1-p,1);
        
        % save results
        save([save_path expt_name{n} '_' expt_ee{e} '_cc_rand_graph_comm_k_' ...
            num2str(k) '.mat'],'num_rand_clique','num_rand_comm',...
            'sz_rand_comm','comm_sz_thresh','-v7.3');
        
    end
end