function [] = fig_core_mc_pot(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
k = param.k;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.ens;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
OSI_thresh = param.OSI_thresh;
num_expt = length(expt_name);
linew = param.linew;
mc_minsz = param.mc_minsz;

qnoise = 0.4;

load(ccode_path);
load(rwbmap);

%% initialize
pred_mLL = cell(num_expt,1);
pred_stats = zeros(2,num_expt,3);
thr = zeros(num_expt,1);

% cos_sim = {};
% cos_sim_avg = [];
% cos_thresh = [];
% pred = {};
% pred_stats = [];

%%
expt_count = 0;
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1)-num_stim;
    num_frame = length(Pks_Frame);
    num_frame_full = size(Spikes,2);
    vis_stim_high = vis_stim(Pks_Frame);
    data_high = Spikes(:,Pks_Frame);
    
    [core_crf,mc_in_core] = find_core_max_clique(best_model.graph,num_stim,mc_minsz);
    
%     % plot model with secondary connections
%     plotHighlightSecondOrder(best_model.graph,Coord_active,num_stim);
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} 'second_order_ensemble.pdf'])

    % prediction
    LL_frame = zeros(num_frame,num_stim);
    for ii = 1:num_frame
        for jj = 1:num_stim
            for kk = 1:length(mc_in_core{jj})
                mc = mc_in_core{jj}{kk}';
                stim_vec = zeros(num_stim,1);
                stim_vec(jj) = 1;
                data_stim = [data_high(mc,ii);stim_vec]';
                mc_ll = compute_avg_log_likelihood(best_model.node_pot...
                    ([mc,num_node+1:num_node+num_stim]),...
                    best_model.edge_pot([mc,num_node+1:num_node+num_stim],...
                    [mc,num_node+1:num_node+num_stim]),best_model.logZ,data_stim);
                LL_frame(ii,jj) = LL_frame(ii,jj)+mc_ll;
            end
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
        true_label = vis_stim_high==ii;
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
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
%         expt_ee '_' ge_type '_all_LL.pdf'])
    
end

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

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats(1,:,1),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,1),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w');
ylabel('Accuracy')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats(1,:,2),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,2),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w');
ylabel('Precision')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats(1,:,3),'positions',1,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
h = boxplot(pred_stats(2,:,3),'positions',1.5,'width',ww,'colors',mycc.blue);
setBoxStyle(h,linew)
xlim([0.5 2]); ylim([0 1])
set(gca,'xcolor','w')
ylabel('Recall')
box off

% print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee '_' ge_type ...
%     '_all_pred_stats.pdf'])


end

