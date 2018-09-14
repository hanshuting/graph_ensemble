% script for figure 7, comparing single model and combined model

%% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high','on_high'}};
num_shuff = 100;
k = 3;
p = 0.05;

% clustering coefficient bin  range
lcc_bin_range = 0:0.02:1;

%% load all data
num_expt = length(expt_name);
comm_sz_all = cell(num_expt,3);
comm_deg_all = cell(num_expt,3);
comm_ov_all = cell(num_expt,3);
comm_mem_all = cell(num_expt,3);
node_deg_all = cell(num_expt,3);
clus_coeff_cum = cell(num_expt,3);

for n = 1:length(expt_name)
    
    % set path
    model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
    cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
    comm_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\comm\'];
    save_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
    fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    model_graph = zeros(num_node,num_node,3);
    model_ep = zeros(num_node,num_node,3);
    
    for e = 1:length(expt_ee)
        
        % load single models
        load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([comm_path expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
            num2str(k) '.mat']);
        num_comm = length(comm_loopy);
        model_graph(:,:,e) = graph;
        model_ep(:,:,e) = edge_pot;
        
        % community size
        comm_sz_all{n,e} = cellfun('length',comm_loopy);
        
        % community degree
        cr_comm_deg = zeros(num_comm,1);
        for i = 1:num_comm
            cr_comm = comm_loopy{i};
            cr_comm_deg(i) = sum(sum(graph(cr_comm,cr_comm)))/2;
        end
        comm_deg_all{n,e} = cr_comm_deg;
        
        % community overlap
        cr_comm_ov = zeros(num_comm,num_comm);
        for i = 1:num_comm
            comm1 = comm_loopy{i};
            for j = 1:num_comm
                comm2 = comm_loopy{j};
                cr_comm_ov(i,j) = length(intersect(comm1,comm2));
            end
        end
        idx = eye(size(cr_comm_ov));
        cr_comm_ov = cr_comm_ov(~idx);
        comm_ov_all{n,e} = cr_comm_ov;
        
        % community membership
        comm_mem_all{n,e} = histc(cell2mat(comm_loopy),1:num_node);
        
        % node degree
        node_deg_all{n,e} = sum(graph,2)/2;
        
        % local clustering coefficient
        lcc = local_cluster_coeff(graph);
        lcc_cum = hist(lcc,lcc_bin_range);
        lcc_cum = lcc_cum/sum(lcc_cum);
        lcc_cum = cumsum(lcc_cum,'reverse');
        clus_coeff_cum{n,e} = lcc_cum;
        
    end
    
    % plot models
    hf = figure;set(hf,'color','w','position',[2160,745,1390,250]);
    set(hf,'PaperSize',[15 3],'PaperPosition',[0 0 15 3]);
    for e = 1:length(expt_ee)
        subplot(1,3,e);
        plotGraphModel(model_graph(:,:,e),Coord_active,model_ep(:,:,e),[]);
        if e==1
            title('horizontal');
        elseif e==2
            title('vertical');
        else
            title('combined');
        end
    end
    
    saveas(hf,[fig_path expt_name{n} '_all_single_combined_model.fig']);
    saveas(hf,[fig_path expt_name{n} '_all_single_combined_model.pdf']);
    
end

%% graph properties
% ------------ comm sz --------------- %
comm_sz_bin_range = 1:max(cell2mat(comm_sz_all(:)));
comm_sz_cum = cellfun(@(x) hist(x,comm_sz_bin_range),comm_sz_all,...
    'uniformoutput',false);
comm_sz_cum = cellfun(@(x) x/sum(x),comm_sz_cum,'uniformoutput',false);
comm_sz_cum = cellfun(@(x) cumsum(x,'reverse'),comm_sz_cum,...
    'uniformoutput',false);

%% plot


