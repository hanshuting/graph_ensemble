function [] = fig2_crf_LL_pred_add_neuron(param)
% THIS VERSION ONLY WORKS FOR 2 STIMULI

% parameters
expt_name = param.expt_name;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.ens;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
num_expt = length(expt_name);
linew = param.linew;

qnoise = 0.3;

load(ccode_path);
load(rwbmap);

%% initialize
pred_mLL = cell(num_expt,1);
pred_stats = zeros(2,num_expt,3);
thr = zeros(num_expt,1);
scores = cell(num_expt,2);
label = cell(num_expt,2);

% cos_sim = {};
% cos_sim_avg = [];
% cos_thresh = [];
% pred = {};
% pred_stats = [];

%% process each experiment
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(setdiff(vis_stim,0));
    num_node = size(best_model.graph,1)-num_stim;
    num_frame = length(Pks_Frame);
    num_frame_full = size(Spikes,2);
    vis_stim_high = vis_stim(Pks_Frame);
    data_high = Spikes(:,Pks_Frame);
    
    % plot model with secondary connections
    plotHighlightSecondOrder(best_model.graph,Coord_active,num_stim);
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_second_order_ensemble.pdf'])
    
    % find second order ensembles
    ens = cell(num_stim,1);
    for ii = 1:num_stim
        fo = setdiff(find(best_model.graph(num_node+ii,:)),num_node+ii:num_node+num_stim)';
        so = cell(length(fo),1);
        for jj = 1:length(fo)
            so{jj} = setdiff(find(best_model.graph(fo(jj),:)),[fo;(num_node+ii:num_node+num_stim)'])';
        end
        ens{ii} = [unique([reshape(unique(cell2mat(so)),[],1);fo]);num_node+ii];
    end
    
%     % prediction with cosine similarity
%     pred_mat = [];
%     for ii = 1:num_stim
%         expt_count = expt_count+1;
%         [pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
%             cos_sim_avg(expt_count,:),acc,prc,rec] = core_cos_sim(setdiff(ens{ii},...
%             num_node+1:num_node+num_stim),data_high',vis_stim_high'==ii);
%         pred_stats(expt_count,:) = [acc,prc,rec];
%         pred_mat(ii,:) = pred{expt_count};
%     end
%     
%     plot_pred(pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
%         vis_stim_high,cmap);
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
%         expt_ee '_' ge_type '_pred_sim.pdf'])
%     
%     % plot prediction
%     plot_pred_raster(pred_mat,vis_stim_high,cmap);
%     print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
%         expt_ee '_' ge_type '_all_pred_raster.pdf'])
    
    % calculate likelihood difference
    LL_frame = zeros(num_frame,num_stim);
    for ii = 1:num_frame
        for jj = 1:num_stim
            stim_vec = zeros(num_stim,1);
            stim_vec(jj) = 1;
            data_stim = [data_high(:,ii);stim_vec]';
            LL_frame(ii,jj) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,data_stim);
        end
    end
    LLs = LL_frame(:,1)-LL_frame(:,2); % currently only works for 2 stim
    
    % threshold by 3 std of noise
    th1 = quantile(LLs(:),qnoise);
    th2 = quantile(LLs(:),1-qnoise);
    LLs_th = LLs;
    LLs_th(LLs<=th1 | LLs>=th2) = NaN;
    thr(n,1) = 3*nanstd(LLs_th(:))+nanmean(LLs_th(:));
    thr(n,2) = -3*nanstd(LLs_th(:))+nanmean(LLs_th(:));
%     sigind = abs(LLs)<=quantile(abs(LLs),qnoise);
%     thr(n) = 3*std(abs(LLs(sigind)))+mean(abs(LLs(sigind)));
    [~,pred] = max(LL_frame,[],2);
    pred(pred==1 & LLs<thr(n,1)) = 0;
    pred(pred==2 & LLs>thr(n,2)) = 0;
    
    % collect LL per stim
    pred_mLL{n,1} = LLs(vis_stim_high==1);
    pred_mLL{n,2} = LLs(vis_stim_high==2);
    
    % prediction statistics
    pred_mat = zeros(num_stim,num_frame);
    pred_mat_full = zeros(num_stim,num_frame_full);
    for ii = 1:num_stim
        pred = reshape(pred,[],1);
        true_label = reshape(vis_stim_high==ii,[],1);
        TP = sum(pred==ii & true_label==1);
        TN = sum(pred~=ii & true_label==0);
        FP = sum(pred==ii & true_label==0);
        FN = sum(pred~=ii & true_label==1);
        acc = (TP+TN)/(TP+TN+FN+FP);
        prc = TP/(TP+FP);
        rec = TP/(TP+FN);
        pred_stats(ii,n,:) = [acc,prc,rec];
        pred_mat(ii,:) = pred==ii;
        pred_mat_full(ii,Pks_Frame) = pred==ii;
    end
    
    % collect scores for ROC
    scores{n,1} = LLs;
    scores{n,2} = -LLs;
    label{n,1} = reshape(vis_stim_high==1,[],1);
    label{n,2} = reshape(vis_stim_high==2,[],1);
    
    % plot prediction
    plot_pred_raster(pred_mat_full,vis_stim,cmap);
%     plot_pred_raster(pred_mat,vis_stim_high,cmap);
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
        expt_ee '_' ge_type '_all_pred_raster.pdf'])
    
    % plot LL
    figure; set(gcf,'color','w','position',[2021 700 807 222])
    set(gcf,'paperpositionmode','auto')
    hold on
    numrep = ceil(max(LLs)-min(LLs));
    imagesc(repmat(vis_stim_high',numrep,1))
    colormap(cmap)
    patch([0 num_frame num_frame 0 0],(numrep+1)/2-[thr(n,1) thr(n,1) thr(n,2)...
        thr(n,2) thr(n,1)],mycc.gray_light,'edgecolor','none','facealpha',0.5);
    plot([0 num_frame],(numrep+1)/2*[1 1],'--','color',mycc.black,'linewidth',linew)
    plot((numrep+1)/2-LLs,'k','linewidth',linew)
    xlim([1 num_frame]); ylim([0.5 0.5+numrep])
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
        expt_ee '_' ge_type '_all_LL.pdf'])
    
end

%% ROC curve
figure;set(gcf,'color','w','position',[2801 517 248 225])
hold on

% individual curve
auc = zeros(num_expt,2);
for ii = 1:num_expt
    auc(ii,1) = plotROCmultic(label{ii,1},scores{ii,1},1,mycc.red_light,linew);
    auc(ii,2) = plotROCmultic(label{ii,2},scores{ii,2},1,mycc.blue_light,linew);
end

% mean curve
plotROCmultic(cell2mat(label(:,1)),cell2mat(scores(:,1)),1,mycc.red,3*linew);
plotROCmultic(cell2mat(label(:,2)),cell2mat(scores(:,2)),1,mycc.blue,3*linew);

set(gca,'linewidth',linew)
legend('off')

saveas(gcf,[fig_path expt_ee '_' ge_type '_roc.pdf'])

%% plot stats
figure;
set(gcf,'color','w','position',[2041 533 993 235]);
set(gcf,'paperpositionmode','auto')

ww = 0.2;

% mean sim value
subplot(1,4,1); hold on
patch([0.5 2 2 0.5 0.5],[thr(1,1) thr(1,1) thr(1,2) thr(1,2) thr(1,1)],...
    mycc.gray_light,'facealpha',0.5,'edgecolor','none')
plot([0.5 2],[0 0],'k--','linewidth',linew)
h = boxplot(pred_mLL{1,1},'positions',1,'width',ww,'colors',mycc.red);
set(h,'linewidth',linew);set(h(7,:),'visible','off');set(h(6,:),'linewidth',2*linew);
h = boxplot(pred_mLL{1,2},'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([-1.5 1.5])
set(gca,'xcolor','w');
ylabel('LL')
box off
pval = ranksum(pred_mLL{1,1},pred_mLL{1,2});
title(num2str(pval));

% AUC
subplot(1,5,2); hold on
h = boxplot(auc(:,1),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(auc(:,2),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
ylabel('AUC')
set(gca,'xcolor','w')
box off
pval = ranksum(pred_stats(1,:,1),pred_stats(2,:,1));
title(num2str(pval));

% accuracy
subplot(1,5,3); hold on
h = boxplot(pred_stats(1,:,1),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,1),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w');
ylabel('Accuracy')
box off
pval = ranksum(pred_stats(1,:,1),pred_stats(2,:,1));
title(num2str(pval));

% precision
subplot(1,5,4); hold on
h = boxplot(pred_stats(1,:,2),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,2),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w');
ylabel('Precision')
box off
pval = ranksum(pred_stats(1,:,3),pred_stats(2,:,3));
title(num2str(pval));

% recall
subplot(1,5,5); hold on
h = boxplot(pred_stats(1,:,3),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,3),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w')
ylabel('Recall')
box off
pval = ranksum(pred_stats(1,:,3),pred_stats(2,:,3));
title(num2str(pval));

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee '_' ge_type ...
    '_all_pred_stats.pdf'])


end

