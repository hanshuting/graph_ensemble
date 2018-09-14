% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis','m37_d2'};
combined_expt = {'on_high_add_neuron','vis_02_high_add_neuron'};
ee = {{'01_high','02_high'},{'vis_02_high'}};
num_shuff = 100;
k = 3;
p = 0.05;


%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm_on\'];
        save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
        % load data
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([cc_path  expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
            num2str(k) '.mat']); % clique_loopy, comm_loopy
        load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_comm_k_' ...
            num2str(k) '.mat']); % clique_cc, comm_cc
        load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_rand_graph_comm_k_' ...
            num2str(k) '.mat']); % clique_cc_rand, comm_cc_rand
        load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
            num2str(k) '.mat']); % clique_loopy_rand, comm_loopy_rand
        
        %% find ensembles
        % ------------- SVD result --------------- %
        ee_core_svd = [];
        for i = 1:length(svd_state)
            if any(expt_ee{e}==svd_state{i})
                ee_core_svd = core_svd{i};
                break;
            end
        end
        
        % ------------- cc graph --------------- %
        % 1. community
        % threshold communities size
        comm_cc_randl = vertcat(comm_cc_rand{:});
        num_comm_cc_rand = length(comm_cc_randl);
        if ~isempty(comm_cc_randl)
            comm_sz_cc_rand = cellfun('length',comm_cc_randl);
            bin_range = 1:max(comm_sz_cc_rand);
            hist_comm_sz_cc_rand = hist(comm_sz_cc_rand,bin_range);
            hist_comm_sz_cc_rand = hist_comm_sz_cc_rand/sum(hist_comm_sz_cc_rand);
            cdf_comm_sz_cc_rand = cumsum(hist_comm_sz_cc_rand);
            comm_sz_thresh_cc = find(cdf_comm_sz_cc_rand>1-p,1);
            comm_sz_cc = cellfun('length',comm_cc);
            comm_cc_thresh = comm_cc(comm_sz_cc>=comm_sz_thresh_cc);
        else
            comm_cc_thresh = comm_cc;
        end
        
        % threshold community membership
        mem_cc = histc(cell2mat(comm_cc_thresh),1:num_node);
        if ~isempty(comm_cc_randl)
            mem_cc_rand = histc(cell2mat(comm_cc_randl),1:num_node)/num_shuff;
            ccdist = fitdist(mem_cc_rand(:),'normal');
            mem_cc_thresh = icdf(ccdist,0.95);
        else
            mem_cc_thresh = 0;
        end
        core_mem_cc = find(mem_cc>mem_cc_thresh);
        core_mem_cc_perc = length(core_mem_cc)/num_node;
        
        % 2. centrality
        cc_cent = eigenvec_centrality(cc_graph);
        cc_rand_cent = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            cc_rand_cent(:,i) = eigenvec_centrality(logical(cc_rand_graph(:,:,i)));
        end
        ccdist = fitdist(cc_rand_cent(:),'normal');
        cc_cent_thresh = icdf(ccdist,1-core_mem_cc_perc);
        core_cent_cc = find(cc_cent>cc_cent_thresh);
        
        % ------------- crf graph --------------- %
        % 1. community
        % threshold communities
        comm_loopy_randl = vertcat(comm_loopy_rand{:});
        num_comm_loopy_rand = length(comm_loopy_randl);
        if ~isempty(comm_loopy_randl)
            comm_sz_loopy_rand = cellfun('length',comm_loopy_randl);
            if length(comm_sz_loopy_rand)<num_shuff
                comm_sz_loopy_rand(end+1:num_shuff) = 0;
            end
            bin_range = 0:max(comm_sz_loopy_rand);
            hist_comm_sz_loopy_rand = hist(comm_sz_loopy_rand,bin_range);
            hist_comm_sz_loopy_rand = hist_comm_sz_loopy_rand/sum(hist_comm_sz_loopy_rand);
            cdf_comm_sz_loopy_rand = cumsum(hist_comm_sz_loopy_rand);
            comm_sz_thresh_loopy = find(cdf_comm_sz_loopy_rand>1-p,1);
            comm_sz_loopy = cellfun('length',comm_loopy);
            comm_loopy_thresh = comm_loopy(comm_sz_loopy>=comm_sz_thresh_loopy);
        else
            comm_loopy_thresh = comm_loopy;
        end
            
        % threshold community membership
        mem_loopy = histc(cell2mat(comm_loopy_thresh),1:num_node);
        if ~isempty(comm_loopy_randl)
            mem_loopy_rand = histc(cell2mat(comm_loopy_randl),1:num_node)/num_shuff;
            ccdist = fitdist(mem_loopy_rand(:),'normal');
            mem_loopy_thresh = icdf(ccdist,0.95);
        else
            mem_loopy_thresh = 0;
        end
        core_mem_loopy = find(mem_loopy>mem_loopy_thresh);
        core_mem_loopy_perc = length(core_mem_loopy)/num_node;

        num_node = size(graph,1);
        num_edge = sum(sum(logical(tril(graph))));
        edge_list = zeros(num_edge,2);
        [edge_list(:,2),edge_list(:,1)] = find(tril(graph));
        G_on = zeros(num_node,num_node);
        for i = 1:size(G,2)
            node_1 = edge_list(i,1);
            node_2 = edge_list(i,2);
            G_on(node_1,node_2) = G(4,i);%-G(3,i)-G(2,i);
            G_on(node_2,node_1) = G(4,i);%-G(3,i)-G(2,i);
        end
        G_on = G_on>0;
%         comm_loopy_on = k_clique(k,G_on);
%         mem_loopy_on = histc(cell2mat(comm_loopy_on),1:num_node);
%         core_mem_loopy = find(mem_loopy_on>0);
%         core_mem_loopy_perc = length(core_mem_loopy)/num_node;
        
        % 2. centrality
        loopy_cent = eigenvec_centrality(G_on);
%         loopy_rand_cent = zeros(num_node,num_shuff);
%         for i = 1:num_shuff
%             loopy_rand_cent(:,i) = eigenvec_centrality(logical(loopy_rand_graph(:,:,i)));
%         end
%         ccdist = fitdist(loopy_rand_cent(:),'normal');
%         loopy_cent_thresh = icdf(ccdist,1-core_mem_loopy_perc);
        loopy_cent_thresh = multithresh(loopy_cent);
        core_cent_loopy = find(loopy_cent>loopy_cent_thresh);
        
        % 3. crf max edge pot subgraph core
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        core_edge_loopy = high_score_subgraph(graph,G);
        
        % 4. add neuron ensemble
        load(['C:\Shuting\fwMatch\results\' expt_name{n} '\models\' expt_name{n} ...
            '_' combined_expt{n} '_loopy_best_model.mat']);
        core_stim_assc = find(graph(num_node+e,:));
        
        % save results
        save([save_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat'],'ee_core_svd','core_mem_cc','core_mem_loopy',...
            'core_cent_cc','core_cent_loopy','core_edge_loopy','core_stim_assc');


        %% plot results
        hf = figure;set(gcf,'color','w','position',[1943,176,1025,573]);
        subplot(2,3,1);
        plotGraphHighlight(Coord_active,ee_core_svd,'k');
        title('SVD');
        subplot(2,3,4);
        plotGraphHighlight(Coord_active,core_edge_loopy,'r');
        title('CRF coact');
        subplot(2,3,2);
        plotGraphHighlight(Coord_active,core_mem_cc,'b');
        title('CC mem');
        subplot(2,3,5);
        plotGraphHighlight(Coord_active,core_mem_loopy,'r');
        title('CRF mem');
        subplot(2,3,3);
        plotGraphHighlight(Coord_active,core_cent_cc,'b');
        title('CC cent');
        subplot(2,3,6);
        plotGraphHighlight(Coord_active,core_cent_loopy,'r');
        title('CC cent');
        
        figure;set(gcf,'color','w','position',[3000 488 288 254])
        plotGraphHighlight(Coord_active,core_stim_assc,'r');
        title('CRF stim assc')
        
        saveas(hf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core.fig']);
        saveas(hf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core.pdf']);
        
    end
end
