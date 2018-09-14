% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high'}};
% ee = {{'01_high','02_high'},{'vis_01_all','vis_02_all','vis_04_all'}};
num_shuff = 100;
k = 4;
p = 0.01;

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
        comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
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
        for i = 1:length(svd_state)
            if any(expt_ee{e}==svd_state{i})
                ee_core_svd = core_svd{i};
                break;
            end
        end
        
        % ------------- cc graph --------------- %
        % 1. clique centrality
        clique_cc_randl = vertcat(clique_cc_rand{:});
        ccent_cc = histc(cell2mat(clique_cc'),1:num_node);
        if ~isempty(clique_cc_randl)
            ccent_cc_rand = histc(cell2mat(clique_cc_randl'),1:num_node)/num_shuff;
            ccdist = fitdist(ccent_cc_rand(:),'normal');
            ccent_cc_thresh = icdf(ccdist,1-p);
        else
            ccent_cc_thresh = 0;
        end
        core_ccent_cc = find(ccent_cc>ccent_cc_thresh);
        
        % 2. eigenvector centrality
        cc_cent = eigenvec_centrality(cc_graph);
        cc_rand_cent = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            cc_rand_cent(:,i) = eigenvec_centrality(logical(cc_rand_graph(:,:,i)));
        end
        ccdist = fitdist(cc_rand_cent(:),'normal');
%         cc_cent_thresh = icdf(ccdist,1-core_mem_cc_perc);
        cc_cent_thresh = icdf(ccdist,0.7);
        core_cent_cc = find(cc_cent>cc_cent_thresh);
        
        % ------------- crf graph --------------- %
        % 1. clique centrality
        %-- remove weak and negative edges --%
%         edge_thresh = quantile(edge_pot(edge_pot>0),0.3);
        edge_thresh = 0;
        graph = edge_pot>edge_thresh;
        [comm_loopy,clique_loopy] = k_clique(k,graph);
        %-- SH 0414 --%
        clique_loopy_randl = vertcat(clique_loopy_rand{:});
        ccent_loopy = histc(cell2mat(clique_loopy'),1:num_node);
        if ~isempty(clique_loopy_randl)
            ccent_loopy_rand = histc(cell2mat(clique_loopy_randl'),1:num_node)/num_shuff;
            ccdist = fitdist(ccent_loopy_rand(:),'normal');
            ccent_loopy_thresh = icdf(ccdist,1-p);
        else
            ccent_loopy_thresh = 0;
        end
        core_ccent_loopy = find(ccent_loopy>ccent_loopy_thresh);
        
        % 2. eigenvector centrality
        trilgraph = tril(graph);
        num_edge = sum(sum(logical(trilgraph)));
        edge_list = zeros(num_edge,2);
        [edge_list(:,2),edge_list(:,1)] = find(trilgraph);
        G_on = zeros(num_node,num_node);
        for i = 1:size(G,2)
            node_1 = edge_list(i,1);
            node_2 = edge_list(i,2);
            G_on(node_1,node_2) = G(4,i);%-G(3,i)-G(2,i);
        end
        
        loopy_cent = eigenvec_centrality(exp(G_on).*double(G_on~=0));
        loopy_rand_cent = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            loopy_rand_cent(:,i) = eigenvec_centrality(logical(loopy_rand_graph(:,:,i)));
        end
        ccdist = fitdist(loopy_rand_cent(:),'normal');
%         loopy_cent_thresh = icdf(ccdist,1-core_mem_loopy_perc);
        loopy_cent_thresh = icdf(ccdist,0.7);
        core_cent_loopy = find(loopy_cent>loopy_cent_thresh);
        
        
        % crf max edge pot subgraph core
        core_edge_loopy = high_score_subgraph(graph,G);
        
        % save results
%          save([save_path expt_name{n} '_' expt_ee{e} '_edge_subgraph_core.mat'], ...
%             'ee_core_svd','core_edge_loopy','core_cent_cc','core_cent_loopy',...
%             'core_ccent_cc','core_ccent_loopy');

        %% plot results
        hf = figure;set(gcf,'color','w','position',[2000,200,1000,580]);
        subplot(2,3,1);
        plotGraphHighlight(Coord_active,ee_core_svd,'k');
        title('SVD');
        subplot(2,3,2)
        plotGraphHighlight(Coord_active,core_ccent_cc,'b');
        title('cc ccent');
        subplot(2,3,5);
        plotGraphHighlight(Coord_active,core_ccent_loopy,'r');
        title('crf ccent');
        subplot(2,3,3)
        plotGraphHighlight(Coord_active,core_cent_cc,'b');
        title('cc ecent');
        subplot(2,3,6);
        plotGraphHighlight(Coord_active,core_cent_loopy,'r');
        title('crf ecent');
        subplot(2,3,4);
        plotGraphHighlight(Coord_active,core_edge_loopy,'r');
        title('crf edge');       
        
%         saveas(hf,[fig_path expt_name{n} '_' expt_ee{e} '_edge_subgraph_core.fig']);
    end
    
end

%% ensemble predictions
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
%         load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
%             '_core.mat']);
        load([core_path expt_name{n} '_' expt_ee{e} '_edge_subgraph_' ...
            '_core.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\Pks_Frames.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} ...
            '_' expt_ee{e} '.mat']);
        load('C:\Shuting\fwMatch\results\wrbmap.mat');
        
        vis_stim_high = vis_stim(Pks_Frame);
        data_high = Spikes(:,Pks_Frame)';
        
        %% ---------- unweighted cosine distance -------%
        ee_core_svd_vec = zeros(num_node,1);
        ee_core_svd_vec(ee_core_svd) = 1;
        sim_svd = pdist2(data_high,ee_core_svd_vec','cosine');
        sim_svd_thresh = 3*std(sim_svd);
        
        core_edge_loopy_vec = zeros(num_node,1);
        core_edge_loopy_vec(core_edge_loopy) = 1;
        sim_edge_loopy = pdist2(data_high,core_edge_loopy_vec','cosine');
        sim_edge_loopy_thresh = 3*std(sim_edge_loopy);
        
        core_cent_loopy_vec = zeros(num_node,1);
        core_cent_loopy_vec(core_cent_loopy) = 1;
        sim_cent_loopy = pdist2(data_high,core_cent_loopy_vec','cosine');
        sim_cent_loopy_thresh = 3*std(sim_cent_loopy);
        
        core_cent_cc_vec = zeros(num_node,1);
        core_cent_cc_vec(core_cent_cc) = 1;
        sim_cent_cc = pdist2(data_high,core_cent_cc_vec','cosine');
        sim_cent_cc_thresh = 3*std(sim_cent_cc);
        
        core_ccent_loopy_vec = zeros(num_node,1);
        core_ccent_loopy_vec(core_ccent_loopy) = 1;
        sim_ccent_loopy = pdist2(data_high,core_ccent_loopy_vec','cosine');
        sim_ccent_loopy_thresh = 3*std(sim_ccent_loopy);
        
        core_ccent_cc_vec = zeros(num_node,1);
        core_ccent_cc_vec(core_ccent_cc) = 1;
        sim_ccent_cc = pdist2(data_high,core_ccent_cc_vec','cosine');
        sim_ccent_cc_thresh = 3*std(sim_ccent_cc);
        
        % change visual stim here
        ee_stim = 1;
        acr_svd = sum(double((1-sim_svd)>sim_svd_thresh)==double...
            (vis_stim_high==ee_stim))/length(vis_stim_high);
        acr_crf = sum(double((1-sim_edge_loopy)>sim_edge_loopy_thresh)==...
            double(vis_stim_high==ee_stim))/length(vis_stim_high);
        acr_cent_cc= sum(double((1-sim_cent_cc)>sim_cent_cc_thresh)==...
            double(vis_stim_high==ee_stim))/length(vis_stim_high);
        acr_cent_loopy = sum(double((1-sim_cent_loopy)>sim_cent_loopy_thresh)==...
            double(vis_stim_high==ee_stim))/length(vis_stim_high);
        
        %% ----------- weighted cosine distance ----------%
        ee_core_svd_vec = zeros(num_node,1);
        ee_core_svd_vec(ee_core_svd) = 1;
        sim_svd = pdist2(data_high,ee_core_svd_vec','cosine');
        
        trilgraph = tril(graph);
        num_edge = sum(sum(logical(trilgraph)));
        edge_list = zeros(num_edge,2);
        [edge_list(:,2),edge_list(:,1)] = find(trilgraph);
        G_on = zeros(num_node,num_node);
        for i = 1:size(G,2)
            node_1 = edge_list(i,1);
            node_2 = edge_list(i,2);
            G_on(node_1,node_2) = G(4,i);%-G(3,i)-G(2,i);
        end
        core_edge_loopy_vec = zeros(num_node,1);
        core_edge_loopy_vec(core_edge_loopy) = mean(exp(G_on(core_edge_loopy,:)),2);
        sim_edge_loopy = pdist2(data_high,core_edge_loopy_vec','cosine');
        
        %% ----- plot svd and max edge subgraph only ----- %
        figure;set(gcf,'color','w','position',[2000,100,1400,900]);
        % svd
        subplot(6,1,1);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_svd+0.5,'k','linewidth',2)
        plot([0,length(sim_svd)],1.5-[sim_svd_thresh,sim_svd_thresh],'b--');
        colormap(cmap)
        title('svd')
        % max edge pot subgraph
        subplot(6,1,2);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_edge_loopy+0.5,'k','linewidth',2)
        plot([0,length(sim_edge_loopy)],1.5-[sim_edge_loopy_thresh,...
            sim_edge_loopy_thresh],'b--');
        colormap(cmap)
        title('crf')
        % cc ccent
        subplot(6,1,3);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_ccent_cc+0.5,'k','linewidth',2)
        plot([0,length(sim_ccent_cc)],1.5-[sim_ccent_cc_thresh,...
            sim_ccent_cc_thresh],'b--');
        colormap(cmap)
        title('cc ccent')
        % loopy ccent
        subplot(6,1,4);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_ccent_loopy+0.5,'k','linewidth',2)
        plot([0,length(sim_ccent_loopy)],1.5-[sim_ccent_loopy_thresh,...
            sim_ccent_loopy_thresh],'b--');
        colormap(cmap)
        title('crf ccent')
        % cc ecent
        subplot(6,1,5);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_cent_cc+0.5,'k','linewidth',2)
        plot([0,length(sim_cent_cc)],1.5-[sim_cent_cc_thresh,...
            sim_cent_cc_thresh],'b--');
        colormap(cmap)
        title('cc ecent')
        % loopy ecent
        subplot(6,1,6);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_cent_loopy+0.5,'k','linewidth',2)
        plot([0,length(sim_cent_loopy)],1.5-[sim_cent_loopy_thresh,...
            sim_cent_loopy_thresh],'b--');
        colormap(cmap)
        title('crf ecent')
        
%         saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
%             '_core_prediction.fig']);
        
    end
    
end
