% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high'}};
vis_stim_seq = {[1,2],[2]};
colorstr = {[1 0 0],[0 0 1]};
num_shuff = 100;
k = 3;
p = 0.05;
nodesz = 50;
rr = 0.7;

M = 3;
N = 3;

%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    ee_stim = vis_stim_seq{n};
    
    coords = [Coord_active(:,1),Coord_active(:,2)];
    hf = figure;
    set(gcf,'color','w','position',[1943,176,1025,573]);
    set(gcf,'PaperPositionMode','auto')
    
    ints_svd = 1:num_node;
    ints_coact = 1:num_node;
    ints_cc_cent = 1:num_node;
    ints_crf_cent = 1:num_node;
    ints_cc_mem = 1:num_node;
    ints_crf_mem = 1:num_node;
    ints_stim_assc = 1:num_node;
    ints_osi = 1:num_node;
    
    for e = 1:length(expt_ee)
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
        % load data
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
        load([core_path expt_name{n} '_core_OSI_' num2str(ee_stim(e)) '.mat']);
        
        num_model = 8;

        %% plot results
        N = ceil(sqrt(num_model));
        M = ceil(num_model/N);
        
        subplot(M,N,1);hold on;
        indx = ee_core_svd;
        ints_svd = intersect(ints_svd,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_svd,1),coords(ints_svd,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('SVD');
        end
        
        subplot(M,N,2);hold on;
        indx = core_stim_assc;
        ints_stim_assc = intersect(ints_stim_assc,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_stim_assc,1),coords(ints_stim_assc,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CRF stim assc');
        end
        
        subplot(M,N,3);hold on;
        indx = core_edge_loopy;
        ints_coact = intersect(ints_coact,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e})
        if e==length(expt_ee)
            scatter(coords(ints_coact,1),coords(ints_coact,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CRF coact');
        end
        
        subplot(M,N,4);hold on;
        indx = core_mem_cc;
        ints_cc_mem = intersect(ints_cc_mem,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_cc_mem,1),coords(ints_cc_mem,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CC comm');
        end
                
        subplot(M,N,5);hold on;
        indx = core_mem_loopy;
        ints_crf_mem = intersect(ints_crf_mem,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_crf_mem,1),coords(ints_crf_mem,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CRF comm');
        end
        
        subplot(M,N,6);hold on;
        indx = core_cent_cc;
        ints_cc_cent = intersect(ints_cc_cent,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_cc_cent,1),coords(ints_cc_cent,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CC cent');
        end

        subplot(M,N,7);hold on;
        indx = core_cent_loopy;
        ints_crf_cent = intersect(ints_crf_cent,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_crf_cent,1),coords(ints_crf_cent,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('CRF cent');
        end
        
        subplot(M,N,8);hold on;
        indx = core_osi;
        ints_osi = intersect(ints_osi,indx);
        scatter(coords(indx,1),coords(indx,2),nodesz,'filled',...
            'markerfacecolor',colorstr{e});
        if e==length(expt_ee)
            scatter(coords(ints_osi,1),coords(ints_osi,2),nodesz,'filled',...
                'markerfacecolor','k');
            scatter(coords(:,1),coords(:,2),nodesz,'k','linewidth',1.5);
            axis off
            title('high OSI');
        end
        
    end
    
    saveas(hf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core.fig']);
    saveas(hf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core.pdf']);
        
end
