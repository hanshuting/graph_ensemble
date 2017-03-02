function [] = fig2_crf_LL_pred_add_neuron_multistim(param)
% This version should work for multiple stimuli
% doesn't work well with LC's datasets, but works for PA datasets

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

load(ccode_path);
load(rwbmap);

%% initialize
pred_mLL = cell(num_expt,1);
pred_stats = zeros(2,num_expt,3);
thr = zeros(num_expt,1);
scores = cell(num_expt,2);
label = cell(num_expt,2);

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
    vis_stim_label = setdiff(unique(vis_stim),0);
    data_high = Spikes(:,Pks_Frame);
    
    % plot model with secondary connections
%     plotHighlightSecondOrder(best_model.graph,Coord_active,num_stim);
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_second_order_ensemble.pdf'])

    % calculate likelihood difference; last column represents no stim
    LL_frame = zeros(num_frame,num_stim+1);
    for ii = 1:num_frame
        for jj = 1:num_stim
            stim_vec = zeros(num_stim,1);
            stim_vec(jj) = 1;
            data_stim = [data_high(:,ii);stim_vec]';
            LL_frame(ii,jj) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,data_stim);
        end
        stim_vec = zeros(num_stim,1);
        data_stim = [data_high(:,ii);stim_vec]';
        LL_frame(ii,num_stim+1) = compute_avg_log_likelihood(best_model.node_pot,...
            best_model.edge_pot,best_model.logZ,data_stim);
    end
    
    % predict by taking the highest LL
    [LL_max,pred] = max(LL_frame,[],2);
    pred(pred==num_stim+1) = 0;
    
    % prediction statistics
    pred_mat = zeros(num_stim,num_frame);
    pred_mat_full = zeros(num_stim,num_frame_full);
    for ii = 1:num_stim
        
        % collect scores for ROC
        scores{n,ii} = LL_frame(:,ii)-LL_max;
        label{n,ii} = reshape(vis_stim_high==vis_stim_label(ii),[],1);
        
        % calculate statistics
        pred = reshape(pred,[],1);
        true_label = reshape(vis_stim_high==vis_stim_label(ii),[],1);
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
    
    % plot prediction
%     plot_pred_raster(pred_mat_full,vis_stim,cmap);
    plot_pred_raster(pred_mat,vis_stim_high,cmap);
%     print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
%         expt_ee '_' ge_type '_all_pred_raster.pdf'])
    
end

%% ROC curve
figure;set(gcf,'color','w','position',[2801 517 248 225])
hold on

% plot individual curves
auc = zeros(num_expt,2);
xx = cell(num_expt,2);
yy = cell(num_expt,2);
for ii = 1:num_expt
    [auc(ii,1),xx{ii,1},yy{ii,1}] = plotROCmultic(label{ii,1},scores{ii,1},1,mycc.red_light,linew);
    [auc(ii,2),xx{ii,2},yy{ii,2}] = plotROCmultic(label{ii,2},scores{ii,2},1,mycc.blue_light,linew);
end

% calculate mean curve
xvec = 0:0.02:1;
ymat = zeros(num_expt,length(xvec),2);
for ii = 1:num_expt
    [~,uid] = unique(xx{ii,1});
    ymat(ii,:,1) = interp1(xx{ii,1}(uid),yy{ii,1}(uid),xvec);
    [~,uid] = unique(xx{ii,2});
    ymat(ii,:,2) = interp1(xx{ii,2}(uid),yy{ii,2}(uid),xvec);
end

% test auc
pval = ranksum(auc(:,1),auc(:,2));
title(num2str(pval));

% plot mean curve
plot(xvec,squeeze(mean(ymat(:,:,1),1)),'color',mycc.red,'linewidth',2*linew)
plot(xvec,squeeze(mean(ymat(:,:,2),1)),'color',mycc.blue,'linewidth',2*linew)
xlim([0 1]); ylim([0 1])
set(gca,'linewidth',linew)
legend('off')

saveas(gcf,[fig_path expt_ee '_' ge_type '_multistim_roc.pdf'])

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
ylabel('accuracy')
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
ylabel('precision')
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
ylabel('recall')
box off
pval = ranksum(pred_stats(1,:,3),pred_stats(2,:,3));
title(num2str(pval));

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee '_' ge_type ...
    '_multistim_pred_stats.pdf'])


end

