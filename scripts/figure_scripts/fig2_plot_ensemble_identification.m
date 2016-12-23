function [] = fig2_plot_ensemble_identification(param)
% identify ensembles by switching on and off each neuron

% parameters
expt_name = param.expt_name;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.ens;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
num_expt = length(expt_name);
linew = param.linew;

qnoise = 0.7;

load(ccode_path);
load(rwbmap);

%% initialize
cos_sim = struct();
cos_sim_avg = struct();
cos_thresh = struct();
pred = struct();
pred_stats = struct();

%%
expt_count = 0;
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1);
    num_frame = length(Pks_Frame);
    num_frame_full = size(Spikes,2);
    vis_stim_high = vis_stim(Pks_Frame);
    data_high = Spikes(:,Pks_Frame);
    
    best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
    best_model.ep_on = best_model.ep_on - tril(best_model.ep_on);
    
    %% find SVD ensemble
    core_svd = cell(num_stim,1);
    for ii = 1:num_stim
        for jj = 1:length(svd_data.svd_state)
            if strcmp(num2str(ii),svd_data.svd_state{jj})
                core_svd{ii} = svd_data.core_svd{jj};
                break;
            end
        end
    end
    
    %% find ensemble with CRF
    % predict each neuron in turn
    LL_frame = zeros(num_node,num_frame,2);
    for ii = 1:num_node
        for jj = 1:num_frame
            frame_vec = data_high(:,jj)';
            frame_vec(ii) = 0;
            LL_frame(ii,jj,1) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,frame_vec);
            frame_vec(ii) = 1;
            LL_frame(ii,jj,2) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,frame_vec);
        end
    end
    LL_on = squeeze(LL_frame(:,:,2)-LL_frame(:,:,1));
    
%     % AUC - two stim version
%     auc = zeros(num_node,num_stim);
%     true_label = double(vis_stim_high==1)';
%     for jj = 1:num_node
%         [~,~,~,auc(jj,1)] = perfcurve(true_label,LL_frame(jj,:,1)-LL_frame(jj,:,2),1);
%     end
%     true_label = double(vis_stim_high==2)';
%     for jj = 1:num_node
%         [~,~,~,auc(jj,2)] = perfcurve(true_label,LL_frame(jj,:,2)-LL_frame(jj,:,1),1);
%     end
%     
%     % find ensembles
%     perc_thresh = 0.4;
%     core_crf = cell(num_stim,1);
%     core_crf{1} = find((auc(:,1)>auc(:,2))&(auc(:,1)>perc_thresh));
%     core_crf{2} = find((auc(:,2)>auc(:,1))&(auc(:,2)>perc_thresh));

    % threshold and count
    thr = zeros(num_node,1);
    LL_pred = nan(num_node,num_frame);
    for ii = 1:num_node
        [LL_pred(ii,:),thr(ii)] = pred_from_LL(LL_on(ii,:),qnoise);
    end
%     LL_pred_count = zeros(num_node,num_stim+1);
%     for ii = 1:num_stim
%         LL_pred_count(:,ii) = nansum(LL_pred(:,vis_stim_high==ii),2)/sum(vis_stim_high==ii);
%     end
%     LL_pred_count(:,end) = nansum(LL_pred(:,vis_stim_high==0),2)/sum(vis_stim_high==0);
%     [max_perc,node_idt] = max(LL_pred_count,[],2);
%     
%     % find percentage threshold
%     perc_vec = 0:0.05:0.5;
%     acc_vec = zeros(length(perc_vec),num_stim);
%     for ii = 1:num_stim
%         true_label = double(vis_stim_high==ii)';
%         for jj = 1:length(perc_vec)
%             core = find(node_idt.*(max_perc>perc_vec(jj))==ii);
%             [~,~,~,~,acc] = core_cos_sim(core,data_high',true_label);
%             acc_vec(jj,ii) = acc;
%         end
%     end
%     [~,perc_thresh] = max(mean(acc_vec,2));
%     perc_thresh = perc_vec(perc_thresh);
%     
%     core_crf = cell(num_stim+1,1);
%     for ii = 1:num_stim+1
%         core_crf{ii} = find(node_idt.*(max_perc>perc_thresh)==ii);
%     end
%     
%     figure; plot(perc_vec,acc_vec);
    
    % calculate TPR and FPR
    TPR = zeros(num_node,num_stim);
    FPR = zeros(num_node,num_stim);
    for ii = 1:num_stim
        for jj = 1:num_node
            TP = sum(LL_pred(jj,:)==1&vis_stim_high'==ii);
            FP = sum(LL_pred(jj,:)~=1&vis_stim_high'==ii);
            TN = sum(LL_pred(jj,:)~=1&vis_stim_high'~=ii);
            FN = sum(LL_pred(jj,:)==1&vis_stim_high'~=ii);
            TPR(jj,ii) = TP/(TP+FN);
            FPR(jj,ii) = FP/(FP+TN);
        end
    end
    
    % find the best threshold for each stimulus
    th_vec = 0.1:0.1:0.9;
    acc_vec = zeros(length(th_vec),num_stim);
    for ii = 1:num_stim
        true_label = vis_stim_high'==ii;
        for jj = 1:length(th_vec)
            core = find(TPR(:,ii)>th_vec(jj)&FPR(:,ii)<th_vec(jj));
            if ~isempty(core)
                [~,~,~,~,acc] = core_cos_sim(core,data_high',true_label);
            else
                acc = 0;
            end
            acc_vec(jj,ii) = acc;
        end
    end
    th = zeros(num_stim,1);
    for ii = 1:num_stim
        [~,best_indx] = max(acc_vec(:,ii));
        th(ii) = th_vec(best_indx);
    end
    
    % identify ensembles
    core_crf = cell(num_stim,1);
%     tprth = 0.4; fprth = 0.4; % hard definition of threshold
    for ii = 1:num_stim
        core_crf{ii} = find(TPR(:,ii)>th(ii)&FPR(:,ii)<th(ii));
    end
    
    % plot each neuron in ROC space
    nodesz = 15;
    figure; set(gcf,'color','w','position',[1984 327 622 274])
    subplot(1,2,1); hold on
    plot([0 1],[0 1],'k--')
    plot([0 th(1)],th(1)*[1 1],'k--')
    plot(th(1)*[1 1],[th(1) 1],'k--')
    scatter(FPR(:,1),TPR(:,1),nodesz,mycc.gray,'filled')
    scatter(FPR(core_crf{1},1),TPR(core_crf{1},1),nodesz,mycc.red,'filled')
    xlim([0 1]); ylim([0 1])
    xlabel('FPR'); ylabel('TPR');
    subplot(1,2,2); hold on
    plot([0 1],[0 1],'k--')
    plot([0 th(2)],th(2)*[1 1],'k--')
    plot(th(2)*[1 1],[th(2) 1],'k--')
    scatter(FPR(:,2),TPR(:,2),nodesz,mycc.gray,'filled')
    scatter(FPR(core_crf{2},2),TPR(core_crf{2},2),nodesz,mycc.blue,'filled')
    xlim([0 1]); ylim([0 1])
    xlabel('FPR'); ylabel('TPR');
    
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_core_ROCspace.pdf'])
    
    %% plot prediction example
%     [~,indx] = max(LL_pred_count(core_crf{1},1));
    [~,indx] = max(TPR(core_crf{1},1)-FPR(core_crf{1},1));
    indx = core_crf{1}(indx);
    
    figure; set(gcf,'color','w','position',[1988 755 232 209])
    plotGraphHighlight(Coord_active,indx,'r');
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
        expt_ee '_cell_' num2str(indx) '.pdf'])
    
    LL_cell_nor = (LL_on(indx,:)-min(LL_on(indx,:)))/(max(LL_on(indx,:))-min(LL_on(indx,:)));
    thr_cell_nor = (thr(indx)-min(LL_on(indx,:)))/(max(LL_on(indx,:))-min(LL_on(indx,:)));
    plot_pred(LL_pred(indx,:),LL_cell_nor,thr_cell_nor,vis_stim_high,cmap)
    set(gcf,'position',[2055 522 1121 127])
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
        expt_ee '_LL_cell_' num2str(indx) '.pdf'])
    
    %% plot count scatter
%     dsz = 20;
%     stepsz = 0.5;
%     cc_dot = {mycc.red,mycc.red_light;mycc.blue,mycc.blue_light};
%     figure; set(gcf,'color','w','position',[2271 387 250 277])
%     hold on
%     noncore = setdiff(1:num_node,cell2mat(core_crf(1:num_stim)))';
%     plot(stepsz*[1;2]*ones(size(noncore))',LL_pred_count(noncore,1:2)',...
%         'color',mycc.gray_light);
%     scatter(stepsz*ones(size(noncore)),LL_pred_count(noncore,1),...
%         dsz,mycc.gray,'filled')
%     scatter(2*stepsz*ones(size(noncore)),LL_pred_count(noncore,2),...
%         dsz,mycc.gray,'filled')
%     for ii = 1:num_stim
%         plot(stepsz*[1;2]*ones(size(core_crf{ii}))',LL_pred_count(core_crf{ii},1:2)',...
%             'color',cc_dot{ii,2});
%         scatter(stepsz*ones(size(core_crf{ii})),LL_pred_count(core_crf{ii},1),...
%             dsz,cc_dot{ii,1},'filled')
%         scatter(2*stepsz*ones(size(core_crf{ii})),LL_pred_count(core_crf{ii},2),...
%             dsz,cc_dot{ii,1},'filled')
%     end
%     plot([0.3,1.2],perc_thresh*[1,1],'k--');
%     xlim([0.3 1.2])
%     set(gca,'xtick',stepsz*[1,2,3],'xticklabel',{'horizontal','vertical'})
%     ylabel('prediction (%)')
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
%         expt_ee '_pred_count.pdf'])
    
    %% plot ensemble highlight
    figure; set(gcf,'color','w','position',[2162 447 434 267])
    subplot(1,2,1)
    plotGraphHighlight(Coord_active,core_crf{1},mycc.red)
    subplot(1,2,2)
    plotGraphHighlight(Coord_active,core_crf{2},mycc.blue)
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
        expt_ee '_vis_core.pdf'])
    
    %% plot SVD+CRF, calc sim
    rr = 1;
    figure;
    set(gcf,'color','w','position',[2041 430 543 338]);
    set(gcf,'paperpositionmode','auto')
    for ii = 1:num_stim

        % plot
        subplot(1,num_stim,ii);
        plotCoreOverlay(Coord_active,core_crf{ii},core_svd{ii},mycc.orange,...
            mycc.green,rr)
        
        true_label = double(vis_stim_high==ii)';
        
        % make it fair
        if ~isempty(core_svd{ii})
            expt_count = expt_count+1;
            
            crf_svd = intersect(core_crf{ii},core_svd{ii});
            num_cell(expt_count) = size(Spikes,1);
            num_crf(expt_count) = length(core_crf{ii});
            num_svd(expt_count) = length(core_svd{ii});
            num_crf_svd(expt_count) = length(crf_svd);
            
            % SVD
            [pred.svd{expt_count},cos_sim.svd{expt_count},cos_thresh.svd{expt_count},...
                cos_sim_avg.svd{expt_count},acc,prc,rec] = core_cos_sim(core_svd{ii},data_high',true_label);
            pred_stats.svd(expt_count,:) = [acc,prc,rec];
            % CRF
            [pred.crf{expt_count},cos_sim.crf{expt_count},cos_thresh.crf{expt_count},...
                cos_sim_avg.crf{expt_count},acc,prc,rec] = core_cos_sim(core_crf{ii},data_high',true_label);
            pred_stats.crf(expt_count,:) = [acc,prc,rec];
        end
        
    end
    
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat'],...
        'core_crf','core_svd');
    
    %% plot cos sim
%     pind = expt_count-num_stim+1:expt_count;
%     sind = [1,3,2,4]; % sort by vis stim type
%     pred_mat = reshape([cell2mat(pred.crf(pind));cell2mat(pred.svd(pind))]',[],num_stim*2)';
%     plot_pred_raster(pred_mat(sind,:),vis_stim_high,cmap)
    
%     print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
%         expt_ee '_mc_svd_core_pred_raster_' ge_type '.pdf'])
    
end

%% plot stats
% figure;
% set(gcf,'color','w','position',[2041 533 993 235]);
% set(gcf,'paperpositionmode','auto')
% 
% stepsz = 0.5;
% binsz = 0.1;
% ww = 0.2;
% 
% % mean sim value
% subplot(1,4,1); hold on
% % svd model
% mcs = cell2mat(cos_sim_avg.svd');
% scatter((stepsz-binsz)*ones(size(mcs(:,1))),mcs(:,1),30,mycc.green_light,'+','linewidth',linew);
% scatter((stepsz+binsz)*ones(size(mcs(:,2))),mcs(:,2),30,mycc.green,'+','linewidth',linew);
% plot([(stepsz-binsz)*ones(size(mcs(:,1))),(stepsz+binsz)*ones(size(mcs(:,1)))]',...
%     mcs','color',mycc.gray);
% plot([stepsz-binsz*1.5,stepsz-binsz*0.5],nanmean(mcs(:,1))*ones(2,1),'color',...
%     mycc.black,'linewidth',linew);
% plot([stepsz+binsz*0.5,stepsz+binsz*1.5],nanmean(mcs(:,2))*ones(2,1),'color',...
%     mycc.black,'linewidth',linew);
% % crf model
% mcs = cell2mat(cos_sim_avg.crf');
% scatter((2*stepsz-binsz)*ones(size(mcs(:,1))),mcs(:,1),30,mycc.orange_light,'+','linewidth',linew);
% scatter((2*stepsz+binsz)*ones(size(mcs(:,2))),mcs(:,2),30,mycc.orange,'+','linewidth',linew);
% plot([(2*stepsz-binsz)*ones(size(mcs(:,1))),(2*stepsz+binsz)*ones(size(mcs(:,1)))]',...
%     mcs','color',mycc.gray);
% plot([2*stepsz-binsz*1.5,2*stepsz-binsz*0.5],mean(mcs(:,1))*ones(2,1),'color',...
%     mycc.black,'linewidth',linew);
% plot([2*stepsz+binsz*0.5,2*stepsz+binsz*1.5],mean(mcs(:,2))*ones(2,1),'color',...
%     mycc.black,'linewidth',linew);
% xlim([0.2 3*stepsz-0.2])
% % ylim([0 1])
% set(gca,'xtick',[1,2]*stepsz);
% set(gca,'xticklabel',{'SVD','CRF'})
% ylabel('Similarity')
% 
% % accuracy
% subplot(1,4,2); hold on
% h = boxplot(pred_stats.svd(:,1),'positions',stepsz,'width',ww,'colors',mycc.green);
% setBoxStyle(h,linew)
% h = boxplot(pred_stats.crf(:,1),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
% setBoxStyle(h,linew)
% xlim([0 3*stepsz]); ylim([0 1])
% set(gca,'xtick',[1,2]*stepsz);
% set(gca,'xticklabel',{'SVD','CRF'})
% ylabel('Accuracy')
% box off
% 
% % precision
% subplot(1,4,3); hold on
% h = boxplot(pred_stats.svd(:,2),'positions',stepsz,'width',ww,'colors',mycc.green);
% setBoxStyle(h,linew)
% h = boxplot(pred_stats.crf(:,2),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
% setBoxStyle(h,linew)
% xlim([0 3*stepsz]); ylim([0 1])
% set(gca,'xtick',[1,2]*stepsz);
% set(gca,'xticklabel',{'SVD','CRF'})
% ylabel('Precision')
% box off
% 
% % recall
% subplot(1,4,4); hold on
% h = boxplot(pred_stats.svd(:,3),'positions',stepsz,'width',ww,'colors',mycc.green);
% setBoxStyle(h,linew)
% h = boxplot(pred_stats.crf(:,3),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
% setBoxStyle(h,linew)
% xlim([0 3*stepsz]); ylim([0 1])
% set(gca,'xtick',[1,2]*stepsz);
% set(gca,'xticklabel',{'SVD','CRF'})
% ylabel('Recall')
% box off
% 
% print(gcf,'-dpdf','-painters',[fig_path 'core_pred_stats.pdf'])


end

