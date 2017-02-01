
% parameters
rng(1000);
expt_name = {'m21_d2_vis','mf2_v1'};
ee = {{'01','02'},{'vis_01_all','vis_02_all','vis_04_all'}};
num_shuff = 100;
k = 3;
p = 0.05;

% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
        save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
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
        
        % threshold communities
        % threshold - loopy
        comm_loopy_randl = vertcat(comm_loopy_rand{:});
        num_comm_loopy_rand = length(comm_loopy_randl);
        comm_sz_loopy_rand = cellfun('length',comm_loopy_randl);
        bin_range = 0:max(comm_sz_loopy_rand);
        hist_comm_sz_loopy_rand = hist(comm_sz_loopy_rand,bin_range);
        hist_comm_sz_loopy_rand = hist_comm_sz_loopy_rand/sum(hist_comm_sz_loopy_rand);
        cdf_comm_sz_loopy_rand = cumsum(hist_comm_sz_loopy_rand);
        comm_sz_thresh_loopy = find(cdf_comm_sz_loopy_rand>1-p,1);
        comm_sz_loopy = cellfun('length',comm_loopy);
        comm_loopy_thresh = comm_loopy(comm_sz_loopy>comm_sz_thresh_loopy);
        % cc
        comm_cc_randl = vertcat(comm_cc_rand{:});
        num_comm_cc_rand = length(comm_cc_randl);
        comm_sz_cc_rand = cellfun('length',comm_cc_randl);
        bin_range = 1:max(comm_sz_cc_rand);
        hist_comm_sz_cc_rand = hist(comm_sz_cc_rand,bin_range);
        hist_comm_sz_cc_rand = hist_comm_sz_cc_rand/sum(hist_comm_sz_cc_rand);
        cdf_comm_sz_cc_rand = cumsum(hist_comm_sz_cc_rand);
        comm_sz_thresh_cc = find(cdf_comm_sz_cc_rand>1-p,1);
        comm_sz_cc = cellfun('length',comm_cc);
        comm_cc_thresh = comm_cc(comm_sz_cc>comm_sz_thresh_cc);
        
        
        %% plot graphs
        h = figure;set(h,'color','w','position',[600 600 900 300]);
        set(h,'PaperSize',[9 3],'PaperPosition',[0 0 9 3]);
        subplot(1,2,1);
        plotGraphHighlight(Coord_active,cell2mat(comm_loopy_thresh),'r')
        title('loopy model communities');colorbar
        subplot(1,2,2);
        plotGraphHighlight(Coord_active,cell2mat(comm_cc_thresh),'b')
        title('corr model communities');colorbar
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_comm_thresh']);
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_comm_thresh.pdf']);
        
        %% plot threshold communities
        h = figure;set(h,'color','w','position',[600 600 900 300]);
        set(h,'PaperSize',[9 3],'PaperPosition',[0 0 9 3]);
        subplot(1,2,1);
        plotGraphModel(graph,Coord_active,edge_pot)
        title('loopy model')
        subplot(1,2,2);
        plotGraphModel(cc_graph,Coord_active,cc_weight)
        title('corr model')
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_graph_model']);
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_graph_model.pdf']);
        
        num_comm_loopy = length(comm_loopy_thresh);
        num_comm_cc = length(comm_cc_thresh);
        
        %% community size
        comm_sz_bin_range = 1:max([comm_sz_loopy;comm_sz_cc;...
            comm_sz_loopy_rand;comm_sz_cc_rand]);
        comm_sz_cum_loopy = hist(comm_sz_loopy,comm_sz_bin_range);
        comm_sz_cum_loopy = comm_sz_cum_loopy/sum(comm_sz_cum_loopy);
        comm_sz_cum_loopy = cumsum(comm_sz_cum_loopy,'reverse');
        comm_sz_cum_cc = hist(comm_sz_cc,comm_sz_bin_range);
        comm_sz_cum_cc = comm_sz_cum_cc/sum(comm_sz_cum_cc);
        comm_sz_cum_cc = cumsum(comm_sz_cum_cc,'reverse');
        comm_sz_cum_loopy_rand = hist(comm_sz_loopy_rand,comm_sz_bin_range);
        comm_sz_cum_loopy_rand = comm_sz_cum_loopy_rand/sum(comm_sz_cum_loopy_rand);
        comm_sz_cum_loopy_rand = cumsum(comm_sz_cum_loopy_rand,'reverse');
        comm_sz_cum_cc_rand = hist(comm_sz_cc_rand,comm_sz_bin_range);
        comm_sz_cum_cc_rand = comm_sz_cum_cc_rand/sum(comm_sz_cum_cc_rand);
        comm_sz_cum_cc_rand = cumsum(comm_sz_cum_cc_rand,'reverse');
        
        %% community degree
        % real data
        comm_deg_loopy = zeros(num_comm_loopy,1);
        for i = 1:num_comm_loopy
            cr_comm = comm_loopy_thresh{i};
            comm_deg_loopy(i) = sum(sum(graph(cr_comm,cr_comm)))/2;
        end
        comm_deg_cc = zeros(num_comm_cc,1);
        for i = 1:num_comm_cc
            cr_comm = comm_cc_thresh{i};
            comm_deg_cc(i) = sum(sum(cc_graph(cr_comm,cr_comm)))/2;
        end
        
        % random data
        comm_deg_loopy_rand = zeros(num_comm_loopy_rand,1);
        count = 0;
        for i = 1:num_shuff
            cr_graph_comm = comm_loopy_rand{i};
            for j = 1:length(cr_graph_comm)
                count = count+1;
                cr_comm = cr_graph_comm{j};
                comm_deg_loopy_rand(count) = sum(sum(loopy_rand_graph(cr_comm,cr_comm,i)))/2;
            end
        end
        comm_deg_cc_rand = zeros(num_comm_cc_rand,1);
        count = 0;
        for i = 1:num_shuff
            cr_graph_comm = comm_cc_rand{i};
            for j = 1:length(cr_graph_comm)
                count = count+1;
                cr_comm = cr_graph_comm{j};
                comm_deg_cc_rand(count) = sum(sum(cc_rand_graph(cr_comm,cr_comm,i)))/2;
            end
        end
        
        % distribution
        comm_deg_bin_range = 0:max([comm_deg_loopy;comm_deg_cc;comm_deg_loopy_rand;comm_deg_cc_rand]);
        comm_deg_cum_loopy = hist(comm_deg_loopy,comm_deg_bin_range);
        comm_deg_cum_loopy = comm_deg_cum_loopy/sum(comm_deg_cum_loopy);
        comm_deg_cum_loopy = cumsum(comm_deg_cum_loopy,'reverse');
        comm_deg_cum_cc = hist(comm_deg_cc,comm_deg_bin_range);
        comm_deg_cum_cc = comm_deg_cum_cc/sum(comm_deg_cum_cc);
        comm_deg_cum_cc = cumsum(comm_deg_cum_cc,'reverse');
        comm_deg_cum_loopy_rand = hist(comm_deg_loopy_rand,comm_deg_bin_range);
        comm_deg_cum_loopy_rand = comm_deg_cum_loopy_rand/sum(comm_deg_cum_loopy_rand);
        comm_deg_cum_loopy_rand = cumsum(comm_deg_cum_loopy_rand,'reverse');
        comm_deg_cum_cc_rand = hist(comm_deg_cc_rand,comm_deg_bin_range);
        comm_deg_cum_cc_rand = comm_deg_cum_cc_rand/sum(comm_deg_cum_cc_rand);
        comm_deg_cum_cc_rand = cumsum(comm_deg_cum_cc_rand,'reverse');
        
        %% community overlap
        % real data - loopy
        comm_ov_loopy = zeros(num_comm_loopy,num_comm_loopy);
        for i = 1:num_comm_loopy
            comm1 = comm_loopy_thresh{i};
            for j = 1:num_comm_loopy
                comm2 = comm_loopy_thresh{j};
                comm_ov_loopy(i,j) = length(intersect(comm1,comm2));
            end
        end
        idx = eye(size(comm_ov_loopy));
        comm_ov_loopy = comm_ov_loopy(~idx);
        
        % real data - cc
        comm_ov_cc = zeros(num_comm_cc,num_comm_cc);
        for i = 1:num_comm_cc
            comm1 = comm_cc_thresh{i};
            for j = 1:num_comm_cc
                comm2 = comm_cc_thresh{j};
                comm_ov_cc(i,j) = length(intersect(comm1,comm2));
            end
        end
        idx = eye(size(comm_ov_cc));
        comm_ov_cc = comm_ov_cc(~idx);
        
        % random graph - loopy
        comm_ov_loopy_rand = [];
        for m = 1:num_shuff
            cr_graph_comm = comm_loopy_rand{m};
            num_cr_comm = length(cr_graph_comm);
            cr_comm_ov = zeros(num_cr_comm,num_cr_comm);
            for i = 1:num_cr_comm
                comm1 = cr_graph_comm{i};
                for j = 1:num_cr_comm
                    comm2 = cr_graph_comm{j};
                    cr_comm_ov(i,j) = length(intersect(comm1,comm2));
                end
            end
            idx = eye(size(cr_comm_ov));
            cr_comm_ov = cr_comm_ov(~idx);
            comm_ov_loopy_rand(end+1:end+length(cr_comm_ov)) = cr_comm_ov;
        end
        comm_ov_loopy_rand = comm_ov_loopy_rand';
        
        % random graph - cc
        comm_ov_cc_rand = [];
        for m = 1:num_shuff
            cr_graph_comm = comm_cc_rand{m};
            num_cr_comm = length(cr_graph_comm);
            cr_comm_ov = zeros(num_cr_comm,num_cr_comm);
            for i = 1:num_cr_comm
                comm1 = cr_graph_comm{i};
                for j = 1:num_cr_comm
                    comm2 = cr_graph_comm{j};
                    cr_comm_ov(i,j) = length(intersect(comm1,comm2));
                end
            end
            idx = eye(size(cr_comm_ov));
            cr_comm_ov = cr_comm_ov(~idx);
            comm_ov_cc_rand(end+1:end+length(cr_comm_ov)) = cr_comm_ov;
        end
        comm_ov_cc_rand = comm_ov_cc_rand';
        
        % distribution
        comm_ov_bin_range = 0:max([comm_ov_loopy;comm_ov_cc;comm_ov_loopy_rand;comm_ov_cc_rand]);
        comm_ov_cum_loopy = hist(comm_ov_loopy,comm_ov_bin_range);
        comm_ov_cum_loopy = comm_ov_cum_loopy/sum(comm_ov_cum_loopy);
        comm_ov_cum_loopy = cumsum(comm_ov_cum_loopy,'reverse');
        comm_ov_cum_cc = hist(comm_ov_cc,comm_ov_bin_range);
        comm_ov_cum_cc = comm_ov_cum_cc/sum(comm_ov_cum_cc);
        comm_ov_cum_cc = cumsum(comm_ov_cum_cc,'reverse');
        comm_ov_cum_loopy_rand = hist(comm_ov_loopy_rand,comm_ov_bin_range);
        comm_ov_cum_loopy_rand = comm_ov_cum_loopy_rand/sum(comm_ov_cum_loopy_rand);
        comm_ov_cum_loopy_rand = cumsum(comm_ov_cum_loopy_rand,'reverse');
        comm_ov_cum_cc_rand = hist(comm_ov_cc_rand,comm_ov_bin_range);
        comm_ov_cum_cc_rand = comm_ov_cum_cc_rand/sum(comm_ov_cum_cc_rand);
        comm_ov_cum_cc_rand = cumsum(comm_ov_cum_cc_rand,'reverse');
        
        %% community membership
        % real data
        comm_mem_loopy = histc(cell2mat(comm_loopy_thresh),1:num_node);
        comm_mem_cc = histc(cell2mat(comm_cc_thresh),1:num_node);
        
        % random graph
        comm_mem_loopy_rand = zeros(num_node,1);
        for i = 1:num_shuff
            comm_mem_loopy_rand = comm_mem_loopy_rand+...
                histc(cell2mat(comm_loopy_rand{i}),1:num_node);
        end
        comm_mem_cc_rand = zeros(num_node,1);
        for i = 1:num_shuff
            comm_mem_cc_rand = comm_mem_cc_rand+...
                histc(cell2mat(comm_cc_rand{i}),1:num_node);
        end
        
        % distribution
        comm_mem_bin_range = 1:max([comm_mem_loopy;comm_mem_cc]);
        comm_mem_cum_loopy = hist(comm_mem_loopy,comm_mem_bin_range);
        comm_mem_cum_loopy = comm_mem_cum_loopy/sum(comm_mem_cum_loopy);
        comm_mem_cum_loopy = cumsum(comm_mem_cum_loopy,'reverse');
        comm_mem_cum_cc = hist(comm_mem_cc,comm_mem_bin_range);
        comm_mem_cum_cc = comm_mem_cum_cc/sum(comm_mem_cum_cc);
        comm_mem_cum_cc = cumsum(comm_mem_cum_cc,'reverse');
        comm_mem_cum_loopy_rand = hist(comm_mem_loopy_rand/num_shuff,comm_mem_bin_range);
        comm_mem_cum_loopy_rand = comm_mem_cum_loopy_rand/sum(comm_mem_cum_loopy_rand);
        comm_mem_cum_loopy_rand = cumsum(comm_mem_cum_loopy_rand,'reverse');
        comm_mem_cum_cc_rand = hist(comm_mem_cc_rand/num_shuff,comm_mem_bin_range);
        comm_mem_cum_cc_rand = comm_mem_cum_cc_rand/sum(comm_mem_cum_cc_rand);
        comm_mem_cum_cc_rand = cumsum(comm_mem_cum_cc_rand,'reverse');
        
        %% save results
        save([save_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_comm_prop.mat'],...
            'comm_sz_bin_range','comm_sz_cum_cc','comm_sz_cum_loopy',...
            'comm_sz_cum_cc_rand','comm_sz_cum_loopy_rand',...
            'comm_deg_bin_range','comm_deg_cum_cc','comm_deg_cum_loopy',...
            'comm_deg_cum_cc_rand','comm_deg_cum_loopy_rand',...
            'comm_ov_bin_range','comm_ov_cum_cc','comm_ov_cum_loopy',...
            'comm_ov_cum_cc_rand','comm_ov_cum_loopy_rand',...
            'comm_mem_bin_range','comm_mem_cum_cc','comm_mem_cum_loopy',...
            'comm_mem_cum_cc_rand','comm_mem_cum_loopy_rand');
        
        %% --------- plot community statistics ---------- %
        ptsz = 10;
        
        h = figure;set(h,'color','w');
        
        % community size
        subplot(2,2,1);hold on;
        scatter(comm_sz_bin_range,comm_sz_cum_cc,ptsz,'filled','b');
        scatter(comm_sz_bin_range,comm_sz_cum_loopy,ptsz,'filled','r');
        scatter(comm_sz_bin_range,comm_sz_cum_cc_rand,ptsz,'filled','k');
        scatter(comm_sz_bin_range,comm_sz_cum_loopy_rand,ptsz,[0.5,0.5,0.5],'filled');
        plot(comm_sz_bin_range,comm_sz_cum_cc,'b');
        plot(comm_sz_bin_range,comm_sz_cum_loopy,'r');
        plot(comm_sz_bin_range,comm_sz_cum_cc_rand,'k');
        plot(comm_sz_bin_range,comm_sz_cum_loopy_rand,'color',[0.5,0.5,0.5]);
        xlim([comm_sz_bin_range(1) comm_sz_bin_range(end)+1]);ylim([0 1])
        xlabel('community size');ylabel('p');
        legend('cc','crf','cc random','crf random');legend boxoff
        
        % community degree
        subplot(2,2,2);hold on;
        scatter(comm_deg_bin_range,comm_deg_cum_cc,ptsz,'filled','b');
        scatter(comm_deg_bin_range,comm_deg_cum_loopy,ptsz,'filled','r');
        scatter(comm_deg_bin_range,comm_deg_cum_cc_rand,ptsz,'filled','k');
        scatter(comm_deg_bin_range,comm_deg_cum_loopy_rand,ptsz,[0.5,0.5,0.5],'filled');
        plot(comm_deg_bin_range,comm_deg_cum_cc,'b');
        plot(comm_deg_bin_range,comm_deg_cum_loopy,'r');
        plot(comm_deg_bin_range,comm_deg_cum_cc_rand,'k');
        plot(comm_deg_bin_range,comm_deg_cum_loopy_rand,'color',[0.5,0.5,0.5]);
        xlim([comm_deg_bin_range(1) comm_deg_bin_range(end)+1]);ylim([0 1])
        xlabel('community degree');ylabel('p');
        legend('cc','crf','cc random','crf random');legend boxoff
        
        % community overlap
        subplot(2,2,3);hold on;
        scatter(comm_ov_bin_range,comm_ov_cum_cc,ptsz,'filled','b');
        scatter(comm_ov_bin_range,comm_ov_cum_loopy,ptsz,'filled','r');
        scatter(comm_ov_bin_range,comm_ov_cum_cc_rand,ptsz,'filled','k');
        scatter(comm_ov_bin_range,comm_ov_cum_loopy_rand,ptsz,[0.5,0.5,0.5],'filled');
        plot(comm_ov_bin_range,comm_ov_cum_cc,'b');
        plot(comm_ov_bin_range,comm_ov_cum_loopy,'r');
        plot(comm_ov_bin_range,comm_ov_cum_cc_rand,'k');
        plot(comm_ov_bin_range,comm_ov_cum_loopy_rand,'color',[0.5,0.5,0.5]);
        xlim([comm_ov_bin_range(1) comm_ov_bin_range(end)]+1);ylim([0 1])
        xlabel('community overlap');ylabel('p');
        legend('cc','crf','cc random','crf random');legend boxoff
        
        % community membership
        subplot(2,2,4);hold on;
        scatter(comm_mem_bin_range,comm_mem_cum_cc,ptsz,'filled','b');
        scatter(comm_mem_bin_range,comm_mem_cum_loopy,ptsz,'filled','r');
        scatter(comm_mem_bin_range,comm_mem_cum_cc_rand,ptsz,'filled','k');
        scatter(comm_mem_bin_range,comm_mem_cum_loopy_rand,ptsz,[0.5,0.5,0.5],'filled');
        plot(comm_mem_bin_range,comm_mem_cum_cc,'b');
        plot(comm_mem_bin_range,comm_mem_cum_loopy,'r');
        plot(comm_mem_bin_range,comm_mem_cum_cc_rand,'k');
        plot(comm_mem_bin_range,comm_mem_cum_loopy_rand,'color',[0.5,0.5,0.5]);
        xlim([comm_mem_bin_range(1) comm_mem_bin_range(end)+1]);ylim([0 1])
        xlabel('community membership');ylabel('p');
        legend('cc','crf','cc random','crf random');legend boxoff
        
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_comm_prop.fig']);
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_comm_prop.pdf']);
        
        %% ------------- graph properties ---------------- %
        % connection density
        dens_cc = sum(cc_graph(:))/num_node/(num_node-1);
        dens_loopy = sum(graph(:))/num_node/(num_node-1);
        dens_cc_rand = sum(cc_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        dens_loopy_rand = sum(loopy_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        
        % node degree
        node_deg_cc = sum(cc_graph,2)/2;
        node_deg_cc_rand = sum(squeeze(sum(cc_rand_graph,3)),2)/2/num_shuff;
        node_deg_loopy = sum(graph,2)/2;
        node_deg_loopy_rand = sum(squeeze(sum(loopy_rand_graph,3)),2)/2/num_shuff;
        % distribution
        node_deg_bin_range = 1:max([node_deg_cc;node_deg_loopy;...
            node_deg_cc_rand;node_deg_loopy_rand]);
        node_deg_cum_cc = hist(node_deg_cc,node_deg_bin_range);
        node_deg_cum_cc = node_deg_cum_cc/sum(node_deg_cum_cc);
        node_deg_cum_cc = cumsum(node_deg_cum_cc,'reverse');
        node_deg_cum_loopy = hist(node_deg_loopy,node_deg_bin_range);
        node_deg_cum_loopy = node_deg_cum_loopy/sum(node_deg_cum_loopy);
        node_deg_cum_loopy = cumsum(node_deg_cum_loopy,'reverse');
        node_deg_cum_cc_rand = hist(node_deg_cc_rand,node_deg_bin_range);
        node_deg_cum_cc_rand = node_deg_cum_cc_rand/sum(node_deg_cum_cc_rand);
        node_deg_cum_cc_rand = cumsum(node_deg_cum_cc_rand,'reverse');
        node_deg_cum_loopy_rand = hist(node_deg_loopy_rand,node_deg_bin_range);
        node_deg_cum_loopy_rand = node_deg_cum_loopy_rand/sum(node_deg_cum_loopy_rand);
        node_deg_cum_loopy_rand = cumsum(node_deg_cum_loopy_rand,'reverse');
        
        % global clustering coefficient
%         gcc_cc = global_cluster_coeff(cc_graph);
%         gcc_loopy = global_cluster_coeff(graph);
%         % cc - random
%         gcc_cc_rand = zeros(num_shuff,1);
%         for i = 1:num_shuff
%             gcc_cc_rand(i) = global_cluster_coeff(cc_graph_rand(:,:,i));
%         end
%         gcc_cc_rand = mean(gcc_cc_rand);
%         % loopy - random
%         gcc_loopy_rand = zeros(num_shuff,1);
%         for i = 1:num_shuff
%             gcc_loopy_rand(i) = global_cluster_coeff(loopy_graph_rand(:,:,i));
%         end
%         gcc_loopy_rand = mean(gcc_loopy_rand);

        % local clustering coefficient
        lcc_bin_range = 0:0.02:1;
        % cc
        lcc_cc = local_cluster_coeff(cc_graph);
        lcc_cum_cc = hist(lcc_cc,lcc_bin_range);
        lcc_cum_cc = lcc_cum_cc/sum(lcc_cum_cc);
        lcc_cum_cc = cumsum(lcc_cum_cc,'reverse');
        % loopy
        lcc_loopy = local_cluster_coeff(graph);
        lcc_cum_loopy = hist(lcc_loopy,lcc_bin_range);
        lcc_cum_loopy = lcc_cum_loopy/sum(lcc_cum_loopy);
        lcc_cum_loopy = cumsum(lcc_cum_loopy,'reverse');
        % cc - random
        lcc_cc_rand = zeros(num_shuff,num_node);
        for i = 1:num_shuff
            lcc_cc_rand(i,:) = local_cluster_coeff(cc_rand_graph(:,:,i));
        end
        lcc_cum_cc_rand = hist(lcc_cc_rand(:),lcc_bin_range);
        lcc_cum_cc_rand = lcc_cum_cc_rand/sum(lcc_cum_cc_rand);
        lcc_cum_cc_rand = cumsum(lcc_cum_cc_rand,'reverse');
        % loopy - random
        lcc_loopy_rand = zeros(num_shuff,num_node);
        for i = 1:num_shuff
            lcc_loopy_rand(i,:) = local_cluster_coeff(loopy_rand_graph(:,:,i));
        end
        lcc_cum_loopy_rand = hist(lcc_loopy_rand(:),lcc_bin_range);
        lcc_cum_loopy_rand = lcc_cum_loopy_rand/sum(lcc_cum_loopy_rand);
        lcc_cum_loopy_rand = cumsum(lcc_cum_loopy_rand,'reverse');

        %% save results
        save([save_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_graph_prop.mat'],...
            'dens_cc','dens_loopy','dens_cc_rand','dens_loopy_rand',...
            'node_deg_cum_cc','node_deg_cum_loopy',...
            'node_deg_cum_cc_rand','node_deg_cum_loopy_rand',...
            'lcc_cc','lcc_loopy','lcc_cc_rand','lcc_loopy_rand','-v7.3');
        
        %% -------------- plot graph properties ------------ %
        h = figure;set(h,'color','w');
        
        % density
        subplot(2,2,1);hold on;
        bar(1,dens_cc_rand,0.5,'facecolor','k');
        bar(2,dens_cc,0.5,'facecolor','b');
        bar(3,dens_loopy_rand,0.5,'facecolor',[0.5 0.5 0.5]);
        bar(4,dens_loopy,0.5,'facecolor','r');
        xlim([0.5 4.5]);
        set(gca,'xtick',1:4,'xticklabel',{'cc','crf','cc rand','crf rand'})
        ylabel('density');box off
        
        % local clustering coefficient
        subplot(2,2,2);hold on;
        scatter(lcc_bin_range,lcc_cum_cc,ptsz,'filled','b');
        scatter(lcc_bin_range,lcc_cum_loopy,ptsz,'filled','r');
        scatter(lcc_bin_range,lcc_cum_cc_rand,ptsz,'filled','k');
        scatter(lcc_bin_range,lcc_cum_loopy_rand,ptsz,[0.5 0.5 0.5],'filled');
        plot(lcc_bin_range,lcc_cum_cc,'b');
        plot(lcc_bin_range,lcc_cum_loopy,'r');
        plot(lcc_bin_range,lcc_cum_cc_rand,'k');
        plot(lcc_bin_range,lcc_cum_loopy_rand,'color',[0.5 0.5 0.5]);
        xlim([lcc_bin_range(1) lcc_bin_range(end)]);ylim([0 1]);
        xlabel('clustering coeff');ylabel('p');box off
        legend('cc','crf','cc rand','crf rand');legend boxoff
        
        % node degree
        subplot(2,2,3);hold on;
        scatter(node_deg_bin_range,node_deg_cum_cc,ptsz,'filled','b');
        scatter(node_deg_bin_range,node_deg_cum_loopy,ptsz,'filled','r');
        scatter(node_deg_bin_range,node_deg_cum_cc_rand,ptsz,'filled','k');
        scatter(node_deg_bin_range,node_deg_cum_loopy_rand,ptsz,[0.5 0.5 0.5],'filled');
        plot(node_deg_bin_range,node_deg_cum_cc,'b');
        plot(node_deg_bin_range,node_deg_cum_loopy,'r');
        plot(node_deg_bin_range,node_deg_cum_cc_rand,'k');
        plot(node_deg_bin_range,node_deg_cum_loopy_rand,'color',[0.5 0.5 0.5]);
        xlim([node_deg_bin_range(1) node_deg_bin_range(end)+1]);ylim([0 1]);
        xlabel('node degree');ylabel('p');
        legend('cc','crf');legend boxoff
        
        % save figure
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_graph_prop.fig']);
        saveas(h,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) '_graph_prop.pdf']);
        
    end

end