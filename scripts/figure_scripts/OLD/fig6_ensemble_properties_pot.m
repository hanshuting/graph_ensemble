% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high'}};
num_stim = 2;
num_rand = 100;
k = 3;
p = 0.05;
qnoise = 0.2;
save_path = 'C:\Shuting\fwMatch\results\stats\';
all_fig_path = 'C:\Shuting\fwMatch\results\fig\';

%% ensemble predictions
vis_stim_seq = {[1,2]};
stim_cstr = {'r','b'};

% initialize
crf_coact_core_stats = cell(length(expt_name),1);
crf_coact_noncore_stats = cell(length(expt_name),1);
crf_coact_rand_stats = cell(length(expt_name),1);

crf_coact_core_ep_nor = cell(length(expt_name),1);
crf_coact_noncore_ep_nor = cell(length(expt_name),1);
crf_coact_rand_ep_nor = cell(length(expt_name),1);

crf_coact_core_pred = cell(length(expt_name),1);
crf_coact_noncore_pred = cell(length(expt_name),1);
crf_coact_rand_pred = cell(length(expt_name),1);

crf_coact_core_ep_avg = cell(length(expt_name),num_stim,num_stim);
crf_coact_noncore_ep_avg = cell(length(expt_name),num_stim,num_stim);
crf_coact_rand_ep_avg = cell(length(expt_name),num_stim,num_stim);

for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    num_expt = length(expt_ee);
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    ee_stim = vis_stim_seq{n};
    
    crf_coact_rand_pred_single = cell(num_expt,num_rand);
    crf_coact_rand_ep_nor_single = cell(num_expt,num_rand);
    
    for e = 1:num_expt
        
        % set path
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\Pks_Frames.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} ...
            '_' expt_ee{e} '.mat']);
        load('C:\Shuting\fwMatch\results\wrbmap.mat');
        
        vis_stim_high = vis_stim(Pks_Frame);
        data_high = Spikes(:,Pks_Frame)';
        
        %% loopy ep subgraph prediction
        % get on state edge potentials
        trilgraph = tril(graph);
        num_edge = sum(sum(logical(trilgraph)));
        edge_list = zeros(num_edge,2);
        [edge_list(:,2),edge_list(:,1)] = find(trilgraph);
        G_on = zeros(size(graph));
        for i = 1:size(G,2)
            node_1 = edge_list(i,1);
            node_2 = edge_list(i,2);
            G_on(node_1,node_2) = G(4,i);
        end
        
        % ---------- ensemble prediction ------------ %
        crf_coact_vec = false(num_node,1);
        crf_coact_vec(core_edge_loopy) = 1;
        pred_ep_core = zeros(size(data_high,1),1);
        for i = 1:length(pred_ep_core)
            fdata = logical(data_high(i,:).*crf_coact_vec');
            pred_ep_core(i) = sum(F(2,fdata))+sum(sum(G_on(fdata,fdata)));
        end
        pred_ep_core_exp = exp(pred_ep_core);
        crf_coact_core_ep_nor{n}{e} = pred_ep_core_exp/max(pred_ep_core_exp);
        
        ts = crf_coact_core_ep_nor{n}{e};
        ep_pot_pred_thresh = 3*std(ts(ts<quantile(ts,qnoise)))+...
            mean(ts(ts<quantile(ts,qnoise)));
        crf_coact_core_pred{n}{e} = crf_coact_core_ep_nor{n}{e}>ep_pot_pred_thresh;
        
        % prediction statistics
        TP = sum(crf_coact_core_pred{n}{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(crf_coact_core_pred{n}{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(crf_coact_core_pred{n}{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(crf_coact_core_pred{n}{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        crf_coact_core_stats{n}(e,:) = [ee_acc,ee_prc,ee_rec];
        
        % ---------- non ensemble prediction ------------ %
        crf_noncoact_vec = true(num_node,1);
        crf_noncoact_vec(core_edge_loopy) = 0;
        pred_ep_noncore = zeros(size(data_high,1),1);
        for i = 1:length(pred_ep_noncore)
            fdata = logical(data_high(i,:).*crf_noncoact_vec');
            pred_ep_noncore(i) = sum(F(2,fdata))+sum(sum(G_on(fdata,fdata)));
        end
        pred_ep_noncore_exp = exp(pred_ep_noncore);
        crf_coact_noncore_ep_nor{n}{e} = pred_ep_noncore_exp/max(pred_ep_noncore_exp);
        
        ts = crf_coact_noncore_ep_nor{n}{e};
        crf_coact_noncore_pred_thresh = 3*std(ts(ts<quantile(ts,qnoise)))+...
            mean(ts(ts<quantile(ts,qnoise)));
        crf_coact_noncore_pred{n}{e} = crf_coact_noncore_ep_nor{n}{e}>...
            crf_coact_noncore_pred_thresh;
        
        % prediction statistics
        TP = sum(crf_coact_noncore_pred{n}{e}==1&vis_stim_high==ee_stim(e));
        TN = sum(crf_coact_noncore_pred{n}{e}==0&vis_stim_high~=ee_stim(e));
        FP = sum(crf_coact_noncore_pred{n}{e}==1&vis_stim_high~=ee_stim(e));
        FN = sum(crf_coact_noncore_pred{n}{e}==0&vis_stim_high==ee_stim(e));
        ee_acc = (TP+TN)/(TP+TN+FN+FP);
        ee_prc = TP/(TP+FP);
        ee_rec = TP/(TP+FN);
        crf_coact_noncore_stats{n}(e,:) = [ee_acc,ee_prc,ee_rec];
        
        % ---------- random cell prediction ------------ %
        num_coact_core = length(core_edge_loopy);
        
        for i = 1:num_rand
        
            crf_rand_vec = false(num_node,1);            
            crf_rand_vec(randperm(num_node,num_coact_core)) = 1;
            pred_ep_rand = zeros(size(data_high,1),1);
            for j = 1:length(pred_ep_rand)
                fdata = logical(data_high(j,:).*crf_rand_vec');
                pred_ep_rand(j) = sum(F(2,fdata))+sum(sum(G_on(fdata,fdata)));
            end
            pred_ep_rand_exp = exp(pred_ep_rand);
            crf_coact_rand_ep_nor_single{e,i} = pred_ep_rand_exp/max(pred_ep_rand_exp);

            ts = crf_coact_rand_ep_nor_single{e,i};
            crf_coact_rand_pred_thresh = 3*std(ts(ts<quantile(ts,qnoise)))+...
                mean(ts(ts<quantile(ts,qnoise)));
            crf_coact_rand_pred_single{e,i} = crf_coact_rand_ep_nor_single{e,i}...
                >crf_coact_rand_pred_thresh;

            % prediction statistics
            TP = sum(crf_coact_rand_pred_single{e,i}==1&vis_stim_high==ee_stim(e));
            TN = sum(crf_coact_rand_pred_single{e,i}==0&vis_stim_high~=ee_stim(e));
            FP = sum(crf_coact_rand_pred_single{e,i}==1&vis_stim_high~=ee_stim(e));
            FN = sum(crf_coact_rand_pred_single{e,i}==0&vis_stim_high==ee_stim(e));
            ee_acc = (TP+TN)/(TP+TN+FN+FP);
            ee_prc = TP/(TP+FP);
            ee_rec = TP/(TP+FN);
            crf_coact_rand_stats{n}((e-1)*num_rand+i,:) = [ee_acc,ee_prc,ee_rec];
        
        end
        
    end
    
    crf_coact_rand_pred{n} = crf_coact_rand_pred_single;
    crf_coact_rand_ep_nor{n} = crf_coact_rand_ep_nor_single;
    
    %% average prediction per stimuli
    for i = 1:num_stim
        for j = 1:num_stim
            crf_coact_core_ep_avg{n,i,j} = crf_coact_core_ep_nor{n}{i}...
                (vis_stim_high==j);
            crf_coact_noncore_ep_avg{n,i,j} = crf_coact_noncore_ep_nor{n}{i}...
                (vis_stim_high==j);
            for k = 1:num_rand
                crf_coact_rand_ep_avg{n,i,j} = [crf_coact_rand_ep_avg{n,i,j};...
                    crf_coact_rand_ep_nor_single{i,k}(vis_stim_high==j)];
            end
        end
    end
    
    %% --------- plot core prediction example ----------%
    num_frame = length(vis_stim_high);
    vis_stim_mat = zeros(num_expt,length(vis_stim_high));
    for i = 1:num_expt
        vis_stim_mat(i,:) = vis_stim_high'.*double(vis_stim_high'==i);
    end

    hf = figure;rr = 0.3;
    set(hf,'color','w','position',[2000,500,800,270],'PaperPositionMode','auto');
    imagesc(vis_stim_mat);
    colormap(cmap);hold on;
    for i = 1:num_expt

        % plot raster
        plot([find(crf_coact_core_pred{n}{i}),find(crf_coact_core_pred{n}{i})]',...
            [(i-0.2)*ones(sum(crf_coact_core_pred{n}{i}),1),...
            i*ones(sum(crf_coact_core_pred{n}{i}),1)]','color',[0,0,0,rr]);

        % plot continuous prediction
        plot(1:num_frame,i+0.4-crf_coact_core_ep_nor{n}{i}*0.4,'color',0*[1,1,1]);
        
    end
    ylim([0.5 2.5])
    set(gca,'xtick',[],'ytick',[])

%     saveas(hf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.fig']);
%     saveas(hf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.pdf']);


end

%% ------------ plot averaged prediction for each stim --------- %
binsz = 0.2;
stepsz = binsz*2/(length(ee_stim)-1);
stepseq = -binsz:stepsz:binsz;
figure;
set(gcf,'color','w','position',[2000,395,650,306],'PaperPositionMode','auto')
for i = 1:num_stim
    subplot(1,num_stim,i);hold on;
    for j = 1:length(ee_stim)
        h = boxplot(cell2mat(crf_coact_core_ep_avg(:,i,j)),...
            'positions',stepseq(j)+1,'width',stepsz*0.5,'colors',stim_cstr{j});
        set(h(7,:),'visible','off')
        h = boxplot(cell2mat(crf_coact_noncore_ep_avg(:,i,j)),...
            'positions',stepseq(j)+2,'width',stepsz*0.5,'colors',stim_cstr{j});
        set(h(7,:),'visible','off')
        h = boxplot(cell2mat(crf_coact_rand_ep_avg(:,i,j)),...
            'positions',stepseq(j)+3,'width',stepsz*0.5,'colors',stim_cstr{j});
        set(h(7,:),'visible','off')
    end
    xlim([0 4]);ylim([0 1]);box off
    ylabel('normalized potential')
    set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
    set(findobj(gca,'type','line'),'linew',1)
    set(gca,'xtick',1:3,'xticklabel',{'CRF coact','CRF noncoact',...
        'CRF rand'},'XTickLabelRotation',45)
end

% saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_crf_coact'...
%     '_core_prediction_stats.fig']);
% saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_crf_coact'...
%     '_core_prediction_stats.pdf']);

%% plot acc, prc, rec
if 0
svd_statl = vertcat(svd_stats{:});
crf_coact_statl = vertcat(crf_coact_stats{:});
cc_comm_statl = vertcat(cc_comm_stats{:});
crf_comm_statl = vertcat(crf_comm_stats{:});
cc_cent_statl = vertcat(cc_cent_stats{:});
crf_cent_statl = vertcat(crf_cent_stats{:});

hf = figure;
set(gcf,'color','w','position',[2265,439,635,260],'PaperPositionMode','auto');
% set(gcf,'PaperSize',[10 3],'PaperPosition',[0 0 10 3]);
ww = 0.3;
posseq = 1.5:0.5:4;
for i = 1:3
    subplot(1,3,i)
    hold on;
    h = boxplot(svd_statl(:,i),'positions',posseq(1),'width',...
        ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_coact_statl(:,i),'positions',posseq(2),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(cc_comm_statl(:,i),'positions',posseq(3),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_comm_statl(:,i),'positions',posseq(4),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(cc_cent_statl(:,i),'positions',posseq(5),...
        'width',ww,'colors','k');
    set(h(7,:),'visible','off')
    h = boxplot(crf_cent_statl(:,i),'positions',posseq(6),...
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
    set(gca,'xtick',posseq,'xticklabel',{'SVD','CRF coact','CC comm',...
        'CRF comm','CC cent','CRF cent'},'XTickLabelRotation',45)
end
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')


saveas(hf,[all_fig_path 'all_ensemble_pred_acc_prc_rec.fig']);
saveas(hf,[all_fig_path 'all_ensemble_pred_acc_prc_rec.pdf']);

%% save results
save([save_path '_all_pred_stats.mat'],'svd_stats','crf_coact_stats',...
    'cc_comm_stats','crf_comm_stats','cc_cent_stats','crf_cent_stats',...
    'crf_coact_pot_stats','-v7.3');

end
