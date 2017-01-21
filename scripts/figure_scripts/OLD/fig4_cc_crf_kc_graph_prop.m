function [] = fig4_cc_crf_prop(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
k = param.k;
p = param.p;
cc_type = param.cc_type;
comm_type = param.comm_type;
fig_path = param.fig_path.graph_prop;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;

comm_sz_bin_range = param.comm_sz_bin_range;
comm_deg_bin_range = param.comm_deg_bin_range;
comm_ov_bin_range = param.comm_ov_bin_range;
comm_mem_bin_range = param.comm_mem_bin_range;
ndeg_bin_range = param.ndeg_bin_range;
lcc_bin_range = param.lcc_bin_range;
cent_bin_range = param.cent_bin_range;

load(ccode_path);

% initialize results
comm_sz = struct();
comm_ov = struct();
comm_deg = struct();
comm_mem = struct();
comm_num = struct();
dens = struct();
lcc = struct();
ndeg = struct();
cent = struct();

comm_sz_cum = struct();
comm_ov_cum = struct();
comm_deg_cum = struct();
comm_mem_cum = struct();
lcc_cum = struct();
ndeg_cum = struct();
cent_cum = struct();

%% go over experiments
expt_count = 0;
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        expt_count = expt_count+1;
        fprintf('Processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set path
        model_path = [result_path_base expt_name{n} '\models\']; 
        cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\']; 
        comm_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
        
        % load data
        best_model = load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([cc_path  expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
            num2str(k) '.mat']); % clique_loopy, comm_loopy
        load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_comm_k_' ...
            num2str(k) '.mat']); % clique_cc, comm_cc
        load([comm_path  expt_name{n} '_' expt_ee{e} '_cc_rand_graph_comm_k_' ...
            num2str(k) '.mat']); % clique_cc_rand, comm_cc_rand
        load([comm_path  expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_comm_k_' ...
            num2str(k) '.mat']); % clique_loopy_rand, comm_loopy_rand
        
        graph = best_model.graph;
        edge_pot = best_model.edge_pot;
        num_node = size(graph,1);
        
        %% threshold communities by size
        % threshold - loopy
        comm_sz(expt_count).loopy = cellfun('length',comm_loopy);
        [comm_sz_thresh_loopy,comm_num(expt_count).loopy_rand,comm_sz(expt_count).loopy_rand] = ...
            thresh_rand_comm(comm_loopy_rand,p);
        comm_loopy_thresh = comm_loopy(comm_sz(expt_count).loopy>comm_sz_thresh_loopy);
        comm_num(expt_count).loopy = length(comm_loopy_thresh);
        
        % cc
        comm_sz(expt_count).cc = cellfun('length',comm_cc);
        [comm_sz_thresh_cc,comm_num(expt_count).cc_rand,comm_sz(expt_count).cc_rand] = ...
            thresh_rand_comm(comm_cc_rand,p);
        comm_cc_thresh = comm_cc(comm_sz(expt_count).cc>comm_sz_thresh_cc);
        comm_num(expt_count).cc = length(comm_cc_thresh);
        
        %% plot graphs
%         figure; set(gcf,'color','w','position',[600 600 900 300],'paperpositionmode','manual');
%         subplot(1,2,1);
%         plotGraphHighlight(Coord_active,cell2mat(comm_loopy_thresh),'k')
%         title('CRF communities');colorbar
%         subplot(1,2,2);
%         plotGraphHighlight(Coord_active,cell2mat(comm_cc_thresh),'k')
%         title('CC communities');colorbar
%         saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_comm_thresh']);
%         saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_comm_thresh.pdf']);
%         
        %% plot threshold communities
        figure;set(gcf,'color','w','position',[600 600 900 300],'paperpositionmode','auto');
        subplot(1,2,1);
        plotGraphModel(graph,[Coord_active;0,0;0,max(Coord_active(:,1))],edge_pot,[])
        title('CRF'); axis equal
        subplot(1,2,2);
        plotGraphModel(cc_graph,[Coord_active;0,0;0,max(Coord_active(:,1))],cc_weight,[])
        title('CC'); axis equal
        saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_graph_model_' savestr]);
        print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_graph_model_' savestr '.pdf'])
%         saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_cc_loopy_graph_model_' savestr '.pdf']);
        
        %% community size
        % cumulative distribution
        comm_sz_cum(expt_count).loopy = calc_cum_dist(comm_sz(expt_count).loopy/...
            num_node,comm_sz_bin_range);
        comm_sz_cum(expt_count).cc = calc_cum_dist(comm_sz(expt_count).cc/...
            num_node,comm_sz_bin_range);
        comm_sz_cum(expt_count).loopy_rand = calc_cum_dist(comm_sz(expt_count).loopy_rand/...
            num_node,comm_sz_bin_range);
        comm_sz_cum(expt_count).cc_rand = calc_cum_dist(comm_sz(expt_count).cc_rand...
            /num_node,comm_sz_bin_range);
        
        %% community degree
        % real data
        comm_deg(expt_count).loopy = cellfun(@(x) sum(sum(graph(x,x)))/2,comm_loopy_thresh)/num_node;
        comm_deg(expt_count).cc = cellfun(@(x) sum(sum(cc_graph(x,x)))/2,comm_cc_thresh)/num_node;
        
        % random data
        comm_deg(expt_count).loopy_rand = [];
        comm_deg(expt_count).cc_rand = [];
        for i = 1:num_shuff
            if ~isempty(comm_loopy_rand{i})
                comm_deg(expt_count).loopy_rand(end+1:end+length(comm_loopy_rand{i})) = ...
                    cellfun(@(x) sum(sum(loopy_rand_graph(x,x,i)))/2,comm_loopy_rand{i});
            end
            if ~isempty(comm_cc_rand{i})
                comm_deg(expt_count).cc_rand(end+1:end+length(comm_cc_rand{i})) = ...
                    cellfun(@(x) sum(sum(cc_rand_graph(x,x,i)))/2,comm_cc_rand{i});
            end
        end
        comm_deg(expt_count).loopy_rand = comm_deg(expt_count).loopy_rand'/num_node;
        comm_deg(expt_count).cc_rand = comm_deg(expt_count).cc_rand'/num_node;
        
        % distribution
        comm_deg_cum(expt_count).loopy = calc_cum_dist(comm_deg(expt_count).loopy,comm_deg_bin_range);
        comm_deg_cum(expt_count).cc = calc_cum_dist(comm_deg(expt_count).cc,comm_deg_bin_range);
        comm_deg_cum(expt_count).loopy_rand = calc_cum_dist(comm_deg(expt_count).loopy_rand,comm_deg_bin_range);
        comm_deg_cum(expt_count).cc_rand = calc_cum_dist(comm_deg(expt_count).cc_rand,comm_deg_bin_range);
        
        %% community overlap
        % real data - loopy
        comm_ov(expt_count).loopy = zeros(comm_num(expt_count).loopy,comm_num(expt_count).loopy);
        for i = 1:comm_num(expt_count).loopy
            comm_ov(expt_count).loopy(i,:) = cellfun(@(x) length(intersect(x,...
                comm_loopy_thresh{i})),comm_loopy_thresh);
        end
        idx = eye(size(comm_ov(expt_count).loopy));
        comm_ov(expt_count).loopy = comm_ov(expt_count).loopy(~idx)/num_node;
        
        % real data - cc
        comm_ov(expt_count).cc = zeros(comm_num(expt_count).cc,comm_num(expt_count).cc);
        for i = 1:comm_num(expt_count).cc
            comm_ov(expt_count).cc(i,:) = cellfun(@(x) length(intersect(x,...
                comm_cc_thresh{i})),comm_cc_thresh);
        end
        idx = eye(size(comm_ov(expt_count).cc));
        comm_ov(expt_count).cc = comm_ov(expt_count).cc(~idx)/num_node;
        
        % random graph - loopy and cc
        comm_ov(expt_count).loopy_rand = [];
        for m = 1:num_shuff
            cr_graph_comm = comm_loopy_rand{m};
            cr_comm_ov = zeros(length(cr_graph_comm),length(cr_graph_comm));
            for i = 1:length(cr_graph_comm)
                cr_comm_ov(i,:) = cellfun(@(x) length(intersect(x,...
                    cr_graph_comm{i})),cr_graph_comm);
            end
            idx = eye(size(cr_comm_ov));
            cr_comm_ov = cr_comm_ov(~idx);
            comm_ov(expt_count).loopy_rand(end+1:end+length(cr_comm_ov)) = cr_comm_ov;
        end
        comm_ov(expt_count).loopy_rand = comm_ov(expt_count).loopy_rand'/num_node;
        
        % random graph - cc
        comm_ov(expt_count).cc_rand = [];
        for m = 1:num_shuff
            cr_graph_comm = comm_cc_rand{m};
            cr_comm_ov = zeros(length(cr_graph_comm),length(cr_graph_comm));
            for i = 1:length(cr_graph_comm)
                cr_comm_ov(i,:) = cellfun(@(x) length(intersect(x,...
                    cr_graph_comm{i})),cr_graph_comm);
            end
            idx = eye(size(cr_comm_ov));
            cr_comm_ov = cr_comm_ov(~idx);
            comm_ov(expt_count).cc_rand(end+1:end+length(cr_comm_ov)) = cr_comm_ov;
        end
        comm_ov(expt_count).cc_rand = comm_ov(expt_count).cc_rand'/num_node;
        
        % distribution
        comm_ov_cum(expt_count).loopy = calc_cum_dist(comm_ov(expt_count).loopy,comm_ov_bin_range);
        comm_ov_cum(expt_count).cc = calc_cum_dist(comm_ov(expt_count).cc,comm_ov_bin_range);
        comm_ov_cum(expt_count).loopy_rand = calc_cum_dist(comm_ov(expt_count).loopy_rand,comm_ov_bin_range);
        comm_ov_cum(expt_count).cc_rand = calc_cum_dist(comm_ov(expt_count).cc_rand,comm_ov_bin_range);
        
        %% community membership
        % real data
        comm_mem(expt_count).loopy = histc(cell2mat(comm_loopy_thresh),1:num_node)/num_node;
        comm_mem(expt_count).cc = histc(cell2mat(comm_cc_thresh),1:num_node)/num_node;
        
        % random graph
        comm_mem(expt_count).loopy_rand = cellfun(@(x) histc(cell2mat(x),1:num_node),...
            comm_loopy_rand,'uniformoutput',false);
        keep_indx = cellfun(@(x) ~isempty(x),comm_mem(expt_count).loopy_rand);
        comm_mem(expt_count).loopy_rand = mean(cell2mat(comm_mem(expt_count).loopy_rand...
            (keep_indx)'),2)/num_node;
        
        comm_mem(expt_count).cc_rand = cellfun(@(x) histc(cell2mat(x),1:num_node),...
            comm_cc_rand,'uniformoutput',false);
        keep_indx = cellfun(@(x) ~isempty(x),comm_mem(expt_count).cc_rand);
        comm_mem(expt_count).cc_rand = mean(cell2mat(comm_mem(expt_count).cc_rand...
            (keep_indx)'),2)/num_node;
        
        % distribution
        comm_mem_cum(expt_count).loopy = calc_cum_dist(comm_mem(expt_count).loopy,comm_mem_bin_range);
        comm_mem_cum(expt_count).cc = calc_cum_dist(comm_mem(expt_count).cc,comm_mem_bin_range);
        comm_mem_cum(expt_count).loopy_rand = calc_cum_dist(comm_mem(expt_count).loopy_rand,comm_mem_bin_range);
        comm_mem_cum(expt_count).cc_rand = calc_cum_dist(comm_mem(expt_count).cc_rand,comm_mem_bin_range);
               
        %% ------------- graph properties ---------------- %
        % 1. connection density
        dens(expt_count).cc = sum(cc_graph(:))/num_node/(num_node-1);
        dens(expt_count).loopy = sum(graph(:))/num_node/(num_node-1);
        dens(expt_count).cc_rand = sum(cc_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        dens(expt_count).loopy_rand = sum(loopy_rand_graph(:))/num_node/(num_node-1)/num_shuff;
        
        % 2. node degree
        ndeg(expt_count).cc = sum(cc_graph,2)/2/num_node;
        ndeg(expt_count).loopy = sum(graph,2)/2/num_node;
        ndeg(expt_count).cc_rand = zeros(num_node,num_shuff);
        ndeg(expt_count).loopy_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            ndeg(expt_count).cc_rand(:,i) = sum(squeeze(cc_rand_graph(:,:,i)),2)/2/num_node;
            ndeg(expt_count).loopy_rand(:,i) = sum(squeeze(loopy_rand_graph(:,:,i)),2)/2/num_node;
        end
        
        % distribution
        ndeg_cum(expt_count).cc = calc_cum_dist(ndeg(expt_count).cc,ndeg_bin_range);
        ndeg_cum(expt_count).loopy = calc_cum_dist(ndeg(expt_count).loopy,ndeg_bin_range);
        ndeg_cum(expt_count).cc_rand = calc_cum_dist(ndeg(expt_count).cc_rand,ndeg_bin_range);
        ndeg_cum(expt_count).loopy_rand = calc_cum_dist(ndeg(expt_count).loopy_rand,ndeg_bin_range);

        % 3. local clustering coefficient
        lcc(expt_count).cc = local_cluster_coeff(cc_graph);
        lcc(expt_count).loopy = local_cluster_coeff(graph);
        lcc(expt_count).cc_rand = zeros(num_node,num_shuff);
        lcc(expt_count).loopy_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            lcc(expt_count).cc_rand(:,i) = local_cluster_coeff(cc_rand_graph(:,:,i));
            lcc(expt_count).loopy_rand(:,i) = local_cluster_coeff(loopy_rand_graph(:,:,i));
        end
        
        lcc_cum(expt_count).cc = calc_cum_dist(lcc(expt_count).cc,lcc_bin_range);
        lcc_cum(expt_count).loopy = calc_cum_dist(lcc(expt_count).loopy,lcc_bin_range);
        lcc_cum(expt_count).cc_rand = calc_cum_dist(lcc(expt_count).cc_rand,lcc_bin_range);
        lcc_cum(expt_count).loopy_rand = calc_cum_dist(lcc(expt_count).loopy_rand,lcc_bin_range);

        % 4. centrality
        cent(expt_count).loopy = eigenvec_centrality(graph);
        cent(expt_count).cc = eigenvec_centrality(cc_graph);
        cent(expt_count).loopy_rand = zeros(num_node,num_shuff);
        cent(expt_count).cc_rand = zeros(num_node,num_shuff);
        for i = 1:num_shuff
            cent(expt_count).loopy_rand(:,i) = eigenvec_centrality(loopy_rand_graph(:,:,i));
            cent(expt_count).cc_rand(:,i) = eigenvec_centrality(cc_rand_graph(:,:,i));
        end
        
        cent_cum(expt_count).loopy = calc_cum_dist(cent(expt_count).loopy,cent_bin_range);
        cent_cum(expt_count).cc = calc_cum_dist(cent(expt_count).cc,cent_bin_range);
        cent_cum(expt_count).loopy_rand = calc_cum_dist(cent(expt_count).loopy_rand,cent_bin_range);
        cent_cum(expt_count).cc_rand = calc_cum_dist(cent(expt_count).cc_rand,cent_bin_range);
        
    end

end

%% save results
save([save_path 'graph_comm_prop_cc_' cc_type '_comm_' comm_type '_k_' ...
    num2str(k) '_' savestr '.mat'],'expt_name','ee',...
    'comm_sz_bin_range','comm_sz','comm_sz_cum','comm_deg_bin_range',...
    'comm_deg','comm_deg_cum','comm_ov_bin_range','comm_ov','comm_ov_cum',...
    'comm_mem_bin_range','comm_mem','comm_mem_cum','dens','ndeg','ndeg_cum',...
    'lcc','lcc_cum','-v7.3');

%% --------- plot community statistics ---------- %
figure;
set(gcf,'color','w','position',[1987 473 571 547],'PaperPositionMode','auto');

% community size
subplot(2,2,1);
plot_graph_prop_cum_single(comm_sz_cum,mycc,comm_sz_bin_range);
xlabel('s^c^o^m');ylabel('p');
legend off; box on

% community degree
subplot(2,2,2);
plot_graph_prop_cum_single(comm_deg_cum,mycc,comm_deg_bin_range);
xlabel('d^c^o^m');ylabel('p');
legend off; box on

% community overlap
subplot(2,2,3);
plot_graph_prop_cum_single(comm_ov_cum,mycc,comm_ov_bin_range);
xlabel('s^o^v');ylabel('p');
legend off; box on

% community membership
subplot(2,2,4);
plot_graph_prop_cum_single(comm_mem_cum,mycc,comm_mem_bin_range);
xlabel('m');ylabel('p');
box on

saveas(gcf,[fig_path 'comm_prop_cc_' cc_type '_comm_' comm_type '_k_' ...
    num2str(k) '_' savestr '.fig']);
saveas(gcf,[fig_path 'comm_prop_cc_' cc_type '_comm_' comm_type '_k_' ...
    num2str(k) '_' savestr '.pdf']);

%% -------------- plot graph properties ------------ %
linew = 1;

figure;
set(gcf,'color','w','position',[1987 473 571 547],'PaperPositionMode','auto');

% density
boxwd = 0.2;
subplot(2,2,1);hold on;
h = boxplot(cell2mat(getNestedField(dens,'cc')),'positions',0.5,'width',...
    boxwd,'colors',mycc.blue);
set(h(7,:),'visible','off')
set(h,'linewidth',2*linew)
h = boxplot(cell2mat(getNestedField(dens,'loopy')),'positions',1,'width',...
    boxwd,'colors',mycc.orange);
set(h(7,:),'visible','off')
set(h,'linewidth',2*linew)
xlim([0 1.5])
ylabel('Density')
set(gca,'xtick',[0.5 1],'xticklabel',{'CC','CRF'},'linewidth',linew)
% xvec = 0.5:0.5:2;    
% scatter(xvec(1)*ones(expt_count,1),cell2mat(getNestedField(dens,'cc_rand')),[],...
%     mycc.purple,'+','linewidth',linew);
% scatter(xvec(2)*ones(expt_count,1),cell2mat(getNestedField(dens,'cc')),[],...
%     mycc.blue,'+','linewidth',linew);
% scatter(xvec(3)*ones(expt_count,1),cell2mat(getNestedField(dens,'loopy_rand')),[],...
%     mycc.black,'+','linewidth',linew);
% scatter(xvec(4)*ones(expt_count,1),cell2mat(getNestedField(dens,'loopy')),[],...
%     mycc.orange,'+','linewidth',linew);
% plot([xvec(1)-mlinesz xvec(1)+mlinesz],nanmean(cell2mat(getNestedField(dens,'cc_rand')))*...
%     ones(1,2),'color',mycc.purple,'linewidth',linew);
% plot([xvec(2)-mlinesz xvec(2)+mlinesz],nanmean(cell2mat(getNestedField(dens,'cc')))*...
%     ones(1,2),'color',mycc.blue,'linewidth',linew);
% plot([xvec(3)-mlinesz xvec(3)+mlinesz],nanmean(cell2mat(getNestedField(dens,'loopy_rand')))*...
%     ones(1,2),'color',mycc.black,'linewidth',linew);
% plot([xvec(4)-mlinesz xvec(4)+mlinesz],nanmean(cell2mat(getNestedField(dens,'loopy')))*...
%     ones(1,2),'color',mycc.orange,'linewidth',linew);
% xlim([0 2.5]);
% set(gca,'xtick',xvec,'xticklabel',{'CC_r_a_n_d','CC','CRF_r_a_n_d','CRF'},...
%     'XTickLabelRotation',45)
% ylabel('density');
% legend off; box on

% local clustering coefficient
subplot(2,2,2);
plot_graph_prop_cum_single(lcc_cum,mycc,lcc_bin_range);
xlabel('clustering coeff');ylabel('p');
legend off; box on

% node degree
subplot(2,2,3);
plot_graph_prop_cum_single(ndeg_cum,mycc,ndeg_bin_range);
xlabel('node degree');ylabel('p');
legend off; box on

% centrality
subplot(2,2,4)
plot_graph_prop_cum_single(cent_cum,mycc,cent_bin_range);
xlabel('centrality')
box on

% save figure
saveas(gcf,[fig_path 'graph_prop_cc_' cc_type '_comm_' comm_type '_k_' ...
    num2str(k) '_' savestr '.fig']);
saveas(gcf,[fig_path 'graph_prop_cc_' cc_type '_comm_' comm_type '_k_' ...
    num2str(k) '_' savestr '.pdf']);

end
