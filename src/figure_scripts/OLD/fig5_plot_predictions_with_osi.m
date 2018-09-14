% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis','m37_d2'};
ee = {{'01_high','02_high'},{'vis_02_high'}};
num_shuff = 100;
num_model = 8;
k = 3;
p = [0.05,0.01];
qnoise = 0.7;
save_path = 'C:\Shuting\fwMatch\results\stats\';
all_fig_path = 'C:\Shuting\fwMatch\results\fig\stats\';

%% ensemble predictions
vis_stim_seq = {[1,2],[2]};
stim_cstr = {'r','b','k'};
num_ee = sum(cellfun('length',ee));

svd_stats = zeros(num_ee,3);
crf_coact_stats = zeros(num_ee,3);
cc_comm_stats = zeros(num_ee,3);
crf_comm_stats = zeros(num_ee,3);
cc_cent_stats = zeros(num_ee,3);
crf_cent_stats = zeros(num_ee,3);
stim_assc_stats = zeros(num_ee,3);
osi_stats = zeros(num_ee,3);

sim_svd_stim = cell(num_ee,1);
sim_edge_loopy_stim = cell(num_ee,1);
sim_cent_cc_stim = cell(num_ee,1);
sim_cent_loopy_stim = cell(num_ee,1);
sim_mem_cc_stim = cell(num_ee,1);
sim_mem_loopy_stim = cell(num_ee,1);
sim_stim_assc_stim = cell(num_ee,1);
sim_osi_stim = cell(num_ee,1);

core_ov_osi = zeros(num_ee,num_model-1);
core_ov_stim = zeros(num_ee,num_model-1);

ee_count = 0;
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    num_expt = length(expt_ee);
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    % initialize
    sim_svd = cell(num_expt,1);
    sim_edge_loopy = cell(num_expt,1);
    sim_cent_loopy = cell(num_expt,1);
    sim_cent_cc = cell(num_expt,1);
    sim_mem_loopy = cell(num_expt,1);
    sim_mem_cc = cell(num_expt,1);
    sim_stim_assc = cell(num_expt,1);
    sim_osi = cell(num_expt,1);
    
    svd_pred = cell(num_expt,1);
    ep_pred = cell(num_expt,1);
    cent_loopy_pred = cell(num_expt,1);
    cent_cc_pred = cell(num_expt,1);
    mem_loopy_pred = cell(num_expt,1);
    mem_cc_pred = cell(num_expt,1);
    stim_assc_pred = cell(num_expt,1);
    osi_pred = cell(num_expt,1);
    
    sim_svd_thresh = zeros(num_expt,1);
    sim_edge_loopy_thresh = zeros(num_expt,1);
    sim_cent_cc_thresh = zeros(num_expt,1);
    sim_cent_loopy_thresh = zeros(num_expt,1);
    sim_mem_cc_thresh = zeros(num_expt,1);
    sim_mem_loopy_thresh = zeros(num_expt,1);
    sim_stim_assc_thresh = zeros(num_expt,1);
    sim_osi_thresh = zeros(num_expt,1);
    
    ee_acc = cell(num_expt,1);
    ee_prc = cell(num_expt,1);
    ee_rec = cell(num_expt,1);
    
    ee_stim = vis_stim_seq{n};
        
    for e = 1:num_expt
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
        load([core_path expt_name{n} '_core_OSI_' num2str(ee_stim(e)) '.mat']);
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\Pks_Frames.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} ...
            '_' expt_ee{e} '.mat']);
        load('C:\Shuting\fwMatch\results\wrbmap.mat');
        
        vis_stim_high = vis_stim(Pks_Frame);
        data_high = Spikes(:,Pks_Frame)';
        
        ee_count = ee_count+1;
        
        core_stim_assc = core_stim_assc';
        core_edge_loopy = core_edge_loopy';
        
        %% ensemble overlaps
        % with high OSI cells
        core_ov_osi(ee_count,1) = length(intersect(core_osi,core_stim_assc))/...
            length(unique([core_osi;core_stim_assc]));
        core_ov_osi(ee_count,2) = length(intersect(core_osi,ee_core_svd))/...
            length(unique([core_osi;ee_core_svd]));
        core_ov_osi(ee_count,3) = length(intersect(core_osi,core_edge_loopy))/...
            length(unique([core_osi;core_edge_loopy]));
        core_ov_osi(ee_count,4) = length(intersect(core_osi,core_mem_cc))/...
            length(unique([core_osi;core_mem_cc]));
        core_ov_osi(ee_count,5) = length(intersect(core_osi,core_mem_loopy))/...
            length(unique([core_osi;core_mem_loopy]));
        core_ov_osi(ee_count,6) = length(intersect(core_osi,core_cent_cc))/...
            length(unique([core_osi;core_cent_cc]));
        core_ov_osi(ee_count,7) = length(intersect(core_osi,core_cent_loopy))/...
            length(unique([core_osi;core_cent_loopy]));
        
        % with stim assc group
        core_ov_stim(ee_count,1) = length(intersect(core_stim_assc,ee_core_svd))/...
            length(unique([core_stim_assc;ee_core_svd]));
        core_ov_stim(ee_count,2) = length(intersect(core_stim_assc,core_edge_loopy))/...
            length(unique([core_stim_assc;core_edge_loopy]));
        core_ov_stim(ee_count,3) = length(intersect(core_stim_assc,core_mem_cc))/...
            length(unique([core_stim_assc;core_mem_cc]));
        core_ov_stim(ee_count,4) = length(intersect(core_stim_assc,core_mem_loopy))/...
            length(unique([core_stim_assc;core_mem_loopy]));
        core_ov_stim(ee_count,5) = length(intersect(core_stim_assc,core_cent_cc))/...
            length(unique([core_stim_assc;core_cent_cc]));
        core_ov_stim(ee_count,6) = length(intersect(core_stim_assc,core_cent_loopy))/...
            length(unique([core_stim_assc;core_cent_loopy]));
        core_ov_stim(ee_count,7) = length(intersect(core_stim_assc,core_osi))/...
            length(unique([core_stim_assc;core_osi]));
        
        %% ---------- cosine similarity -------%
        ee_core_svd_vec = zeros(num_node,1);
        ee_core_svd_vec(ee_core_svd) = 1;
        sim_svd{e} = 1-pdist2(data_high,ee_core_svd_vec','cosine');
        ts = sim_svd{e};
        sim_svd_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_svd_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_svd_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        svd_pred{e} = ts>sim_svd_thresh(e);
        
        core_edge_loopy_vec = zeros(num_node,1);
        core_edge_loopy_vec(core_edge_loopy) = 1;
        sim_edge_loopy{e} = 1-pdist2(data_high,core_edge_loopy_vec','cosine');
        ts = sim_edge_loopy{e};
        sim_edge_loopy_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_edge_loopy_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_edge_loopy_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        ep_pred{e} = ts>sim_edge_loopy_thresh(e);
        
        core_cent_loopy_vec = zeros(num_node,1);
        core_cent_loopy_vec(core_cent_loopy) = 1;
        sim_cent_loopy{e} = 1-pdist2(data_high,core_cent_loopy_vec','cosine');
        ts = sim_cent_loopy{e};
        sim_cent_loopy_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_cent_loopy_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_cent_loopy_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        cent_loopy_pred{e} = ts>sim_cent_loopy_thresh(e);
        
        core_cent_cc_vec = zeros(num_node,1);
        core_cent_cc_vec(core_cent_cc) = 1;
        sim_cent_cc{e} = 1-pdist2(data_high,core_cent_cc_vec','cosine');
        ts = sim_cent_cc{e};
        sim_cent_cc_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_cent_cc_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_cent_cc_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        cent_cc_pred{e} = ts>sim_cent_cc_thresh(e);
        
        core_mem_loopy_vec = zeros(num_node,1);
        core_mem_loopy_vec(core_mem_loopy) = 1;
        sim_mem_loopy{e} = 1-pdist2(data_high,core_mem_loopy_vec','cosine');
        ts = sim_mem_loopy{e};
        sim_mem_loopy_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_mem_loopy_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_mem_loopy_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        mem_loopy_pred{e} = ts>sim_mem_loopy_thresh(e);
        
        core_mem_cc_vec = zeros(num_node,1);
        core_mem_cc_vec(core_mem_cc) = 1;
        sim_mem_cc{e} = 1-pdist2(data_high,core_mem_cc_vec','cosine');
        ts = sim_mem_cc{e};
        sim_mem_cc_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_mem_cc_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_mem_cc_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        mem_cc_pred{e} = ts>sim_mem_cc_thresh(e);
        
        core_stim_assc_vec = zeros(num_node,1);
        core_stim_assc_vec(core_stim_assc) = 1;
        sim_stim_assc{e} = 1-pdist2(data_high,core_stim_assc_vec','cosine');
        ts = sim_stim_assc{e};
        sim_stim_assc_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_stim_assc_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_stim_assc_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        stim_assc_pred{e} = ts>sim_stim_assc_thresh(e);
        
        core_osi_vec = zeros(num_node,1);
        core_osi_vec(core_osi) = 1;
        sim_osi{e} = 1-pdist2(data_high,core_osi_vec','cosine');
        ts = sim_osi{e};
        sim_osi_stim{ee_count,1} = ts(vis_stim_high~=ee_stim(e));
        sim_osi_stim{ee_count,2} = ts(vis_stim_high==ee_stim(e));
        sim_osi_thresh(e) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
            mean(ts(ts<=quantile(ts,qnoise)));
        osi_pred{e} = ts>sim_osi_thresh(e);
        
        %% cosine similarity stats
        % svd
        TP = sum(svd_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(svd_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(svd_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(svd_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        svd_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % crf coact
        TP = sum(ep_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(ep_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(ep_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(ep_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        crf_coact_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % cc mem
        TP = sum(mem_cc_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(mem_cc_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(mem_cc_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(mem_cc_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        cc_comm_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % crf mem
        TP = sum(mem_loopy_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(mem_loopy_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(mem_loopy_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(mem_loopy_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        crf_comm_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % cc cent
        TP = sum(cent_cc_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(cent_cc_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(cent_cc_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(cent_cc_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        cc_cent_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % crf cent
        TP = sum(cent_loopy_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(cent_loopy_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(cent_loopy_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(cent_loopy_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        crf_cent_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % crf add neuron model
        TP = sum(stim_assc_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(stim_assc_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(stim_assc_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(stim_assc_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        stim_assc_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
        % high OSI cells
        TP = sum(osi_pred{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(osi_pred{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(osi_pred{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(osi_pred{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        osi_stats(ee_count,:) = [ee_acc,ee_prc,ee_rec];
        
    end

    
    %% --------- plot similarity and prediction raster ----------%
    num_frame = length(vis_stim_high);
    num_model = 8;
    vis_stim_mat = repmat(vis_stim_high,1,num_expt*num_model*2)';
    
    figure;rr = 1;ss = 0.1;
    set(gcf,'color','w','position',[1982,15,719,964],'PaperPositionMode','auto');
    imagesc(vis_stim_mat);
    colormap(cmap);hold on;
    for i = 1:num_expt
        
        % plot raster
        plot([find(stim_assc_pred{i}),find(stim_assc_pred{i})]',...
            [(0.5+ss+(i-1)*num_model/2)*ones(sum(stim_assc_pred{i}),1),...
            (1-ss+(i-1)*num_model/2)*ones(sum(stim_assc_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(svd_pred{i}),find(svd_pred{i})]',...
            [(1+ss+(i-1)*num_model/2)*ones(sum(svd_pred{i}),1),...
            (1.5-ss+(i-1)*num_model/2)*ones(sum(svd_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(ep_pred{i}),find(ep_pred{i})]',...
            [(1.5+ss+(i-1)*num_model/2)*ones(sum(ep_pred{i}),1),...
            (2-ss+(i-1)*num_model/2)*ones(sum(ep_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(mem_cc_pred{i}),find(mem_cc_pred{i})]',...
            [(2+ss+(i-1)*num_model/2)*ones(sum(mem_cc_pred{i}),1),...
            (2.5-ss+(i-1)*num_model/2)*ones(sum(mem_cc_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(mem_loopy_pred{i}),find(mem_loopy_pred{i})]',...
            [(2.5+ss+(i-1)*num_model/2)*ones(sum(mem_loopy_pred{i}),1),...
            (3-ss+(i-1)*num_model/2)*ones(sum(mem_loopy_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(cent_cc_pred{i}),find(cent_cc_pred{i})]',...
            [(3+ss+(i-1)*num_model/2)*ones(sum(cent_cc_pred{i}),1),...
            (3.5-ss+(i-1)*num_model/2)*ones(sum(cent_cc_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(cent_loopy_pred{i}),find(cent_loopy_pred{i})]',...
            [(3.5+ss+(i-1)*num_model/2)*ones(sum(cent_loopy_pred{i}),1),...
            (4-ss+(i-1)*num_model/2)*ones(sum(cent_loopy_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot([find(osi_pred{i}),find(osi_pred{i})]',...
            [(4+ss+(i-1)*num_model/2)*ones(sum(osi_pred{i}),1),...
            (4.5-ss+(i-1)*num_model/2)*ones(sum(osi_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        
        % plot sim
        plot(1:num_frame,1.5+(i+num_expt/2-1)*num_model-sim_stim_assc{i},'color',0*[1,1,1])
        plot([1,num_frame],1.5+(i+num_expt/2-1)*num_model-sim_stim_assc_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,2.5+(i+num_expt/2-1)*num_model-sim_svd{i},'color',0*[1,1,1]);
        plot([1,num_frame],2.5+(i+num_expt/2-1)*num_model-sim_svd_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,3.5+(i+num_expt/2-1)*num_model-sim_edge_loopy{i},'color',0*[1,1,1])
        plot([1,num_frame],3.5+(i+num_expt/2-1)*num_model-sim_edge_loopy_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,4.5+(i+num_expt/2-1)*num_model-sim_mem_cc{i},'color',0*[1,1,1])
        plot([1,num_frame],4.5+(i+num_expt/2-1)*num_model-sim_mem_cc_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,5.5+(i+num_expt/2-1)*num_model-sim_mem_loopy{i},'color',0*[1,1,1])
        plot([1,num_frame],5.5+(i+num_expt/2-1)*num_model-sim_mem_loopy_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,6.5+(i+num_expt/2-1)*num_model-sim_cent_cc{i},'color',0*[1,1,1])
        plot([1,num_frame],6.5+(i+num_expt/2-1)*num_model-sim_cent_cc_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,7.5+(i+num_expt/2-1)*num_model-sim_cent_loopy{i},'color',0*[1,1,1])
        plot([1,num_frame],7.5+(i+num_expt/2-1)*num_model-sim_cent_loopy_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        plot(1:num_frame,8.5+(i+num_expt/2-1)*num_model-sim_osi{i},'color',0*[1,1,1])
        plot([1,num_frame],8.5+(i+num_expt/2-1)*num_model-sim_osi_thresh(i)*[1,1],...
            '-','color',0.8*[1,1,1]);
        
    end
    ylim([0.5 num_model*num_expt*1.5+0.5])
    set(gca,'tickdir','out','ytick',[])
    
    saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.fig']);
    saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.pdf']);

end

%% plot average similarity
% significance test
pval = zeros(num_model,1);
[~,pval(1)] = ttest2(cell2mat(sim_stim_assc_stim(:,1)),cell2mat(sim_stim_assc_stim(:,2)));
[~,pval(2)] = ttest2(cell2mat(sim_svd_stim(:,1)),cell2mat(sim_svd_stim(:,2)));
[~,pval(3)] = ttest2(cell2mat(sim_edge_loopy_stim(:,1)),cell2mat(sim_edge_loopy_stim(:,2)));
[~,pval(4)] = ttest2(cell2mat(sim_mem_cc_stim(:,1)),cell2mat(sim_mem_cc_stim(:,2)));
[~,pval(5)] = ttest2(cell2mat(sim_mem_loopy_stim(:,1)),cell2mat(sim_mem_loopy_stim(:,2)));
[~,pval(6)] = ttest2(cell2mat(sim_cent_cc_stim(:,1)),cell2mat(sim_cent_cc_stim(:,2)));
[~,pval(7)] = ttest2(cell2mat(sim_cent_loopy_stim(:,1)),cell2mat(sim_cent_loopy_stim(:,2)));
[~,pval(8)] = ttest2(cell2mat(sim_osi_stim(:,1)),cell2mat(sim_osi_stim(:,2)));

binsz = 0.2;
cstr = {[0 0 0],[1 0.5 0]};
stepsz = binsz*2/(length(ee_stim));
stepseq = -binsz:stepsz:binsz;
figure;
set(gcf,'color','w','position',[2000,441,354,260],'PaperPositionMode','auto')
hold on;
for i = 1:2
    
    indx = 1;
    h = boxplot(cell2mat(sim_stim_assc_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 2;
    h = boxplot(cell2mat(sim_svd_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 3;
    h = boxplot(cell2mat(sim_edge_loopy_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 4;
    h = boxplot(cell2mat(sim_mem_cc_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 5;
    h = boxplot(cell2mat(sim_mem_loopy_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 6;
    h = boxplot(cell2mat(sim_cent_cc_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 7;
    h = boxplot(cell2mat(sim_cent_loopy_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
    indx = 8;
    h = boxplot(cell2mat(sim_osi_stim(:,i)),...
        'positions',stepseq(i)+indx,'width',stepsz*0.5,'colors',cstr{i});
    set(h(7,:),'visible','off')
    if i==length(ee_stim)&&pval(indx)<=p(2)
        scatter(indx+0.01*[1,-1],0.9*[1,1],'k*')
    elseif i==length(ee_stim)&&pval(indx)<=p(1)
        scatter(indx,0.9,'k*')
    end
    
end
xlim([0 num_model+1]);ylim([0 1]);box off
ylabel('similarity')
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'type','line'),'linew',1)
set(gca,'xtick',1:num_model,'xticklabel',{'CRF stim assc','SVD','CRF coact','CC comm',...
    'CRF comm','CC cent','CRF cent','high OSI'},'XTickLabelRotation',45)

saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_sim_avg.fig']);
saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_sim_avg.pdf']);

%% plot acc, prc, rec
hf = figure;
set(gcf,'color','w','position',[2265,439,928,260],'PaperPositionMode','auto');
ww = 0.3;
posseq = 1.5:0.5:(1+num_model*0.5);
for i = 1:3
    subplot(1,3,i)
    hold on;
    h = boxplot(stim_assc_stats(:,i),'positions',posseq(1),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(svd_stats(:,i),'positions',posseq(2),'width',...
        ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_coact_stats(:,i),'positions',posseq(3),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(cc_comm_stats(:,i),'positions',posseq(4),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_comm_stats(:,i),'positions',posseq(5),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(cc_cent_stats(:,i),'positions',posseq(6),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_cent_stats(:,i),'positions',posseq(7),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(osi_stats(:,i),'positions',posseq(8),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    
    xlim([posseq(1)-ww*1.5 posseq(end)+ww*1.5]);ylim([0 1]);box off
    if i==1
        ylabel('accuracy')
    elseif i==2
        ylabel('precision')
    else
        ylabel('recall')
    end
    set(findobj(gca,'type','line'),'linew',1)
    set(gca,'xtick',posseq,'xticklabel',{'CRF stim assc','SVD','CRF coact','CC comm',...
        'CRF comm','CC cent','CRF cent','high OSI'},'XTickLabelRotation',45)
end
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

saveas(hf,[all_fig_path 'all_models_pred_acc_prc_rec.fig']);
saveas(hf,[all_fig_path 'all_models_pred_acc_prc_rec.pdf']);

%% plot core overlaps
hf = figure;
set(gcf,'color','w','position',[2265,430,517,269],'PaperPositionMode','auto');
ww = 0.3;
posseq = 1.5:0.5:(0.5+num_model*0.5);

% overlap with high OSI
subplot(1,2,1)
hold on;
for i = 1:num_model-1
    h = boxplot(core_ov_osi(:,i),'positions',posseq(i),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
end
xlim([posseq(1)-ww*1.5 posseq(end)+ww*1.5]);ylim([0 1]);box off
ylabel('% overlap with high OSI')
set(findobj(gca,'type','line'),'linew',1)
set(gca,'xtick',posseq,'xticklabel',{'CRF stim assc','SVD','CRF coact','CC comm',...
    'CRF comm','CC cent','CRF cent'},'XTickLabelRotation',45)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

% overlap with CRF stim assc
subplot(1,2,2)
hold on;
for i = 1:num_model-1
    h = boxplot(core_ov_stim(:,i),'positions',posseq(i),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
end
xlim([posseq(1)-ww*1.5 posseq(end)+ww*1.5]);ylim([0 1]);box off
ylabel('% overlap with CRF stim assc')
set(findobj(gca,'type','line'),'linew',1)
set(gca,'xtick',posseq,'xticklabel',{'SVD','CRF coact','CC comm',...
    'CRF comm','CC cent','CRF cent','high OSI'},'XTickLabelRotation',45)
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')

saveas(hf,[all_fig_path 'core_overlap.fig']);
saveas(hf,[all_fig_path 'core_overlap.pdf']);

%% save results
save([save_path '_all_pred_stats.mat'],'svd_stats','crf_coact_stats',...
    'cc_comm_stats','crf_comm_stats','cc_cent_stats','crf_cent_stats',...
    'stim_assc_stats','osi_stats','-v7.3');


