function [] = fig6_ensemble_reduction(param)

% parse parameters
expt_name = param.expt_name;
ee = param.ee;
data_path = param.data_path;
fig_path = param.fig_path.stats;
ge_type = param.ge_type;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
linew = param.linew;

bkrmap = load(param.bkrmap);
load(ccode_path);

sample_step = 0.1;
num_rand = 100;

%% ensemble predictions
num_expt = length(expt_name);
sample_seq = -0.9:sample_step:0.9;
[~,indx] = min(abs(sample_seq));
sample_seq(indx) = 0;
rand_seq = 0.1:sample_step:1;

% initialize
stats_pl_all = [];
stats_rd_all = [];

msim_pl_all = [];
msim_rd_all = [];

sim_pl_all = {};
sim_rd_all = {};

true_label_all = {};

%% predict with random cores
expt_count = 0;
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1)-num_stim;
    
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    load([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat'],...
        'core_crf','core_svd');
    
    %% predict with random ensembles
    for s = 1:num_stim
        
        expt_count = expt_count+1;
        ensemble = core_crf{s};
        num_core = length(ensemble);
        noncore = setdiff(1:num_node,ensemble);
        core_plus_seq = round(num_core*sample_seq);
        
        true_label = double(vis_stim_high==s)';
        true_label_all{expt_count} = true_label;
        
        % ensemble plus/minus
        for ii = 1:length(sample_seq)
            for jj = 1:num_rand
                
                % randomly add or remove ensemble neurons
                if core_plus_seq(ii) < 0
                    rand_core = ensemble(randperm(num_core,num_core+core_plus_seq(ii)));
                else
                    rand_core = noncore(randperm(length(noncore),core_plus_seq(ii)));
                    rand_core = [ensemble;rand_core'];
                end
                
                % predict
                [~,sim_core,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data,true_label);
                
                sim_pl_all{expt_count,ii,jj} = sim_core;
                msim_pl_all(expt_count,ii,jj,:) = sim_avg;
                stats_pl_all(expt_count,ii,jj,:) = [acc,prc,rec];
                
            end
        end
    
        % random ensembles
        core_rand_seq = round(num_node*rand_seq);
        for ii = 1:length(rand_seq)
            for jj = 1:num_rand
                
                % randomly sample neurons
                rand_core = randperm(num_node,core_rand_seq(ii));
                
                % predict
                [~,sim_core,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data,true_label);
                
                sim_rd_all{expt_count,ii,jj} = sim_core;
                msim_rd_all(expt_count,ii,jj,:) = sim_avg;
                stats_rd_all(expt_count,ii,jj,:) = [acc,prc,rec];
                
            end
        end
    end
    
end

%% plot 
binsz = 0.02;
wr = 0.3;

figure
set(gcf,'color','w')
set(gcf,'position',[2041 217 649 756]);
set(gcf,'paperpositionmode','auto')

% plus/minus similarity
subplot(4,2,1);
gcapos = get(gca,'position');
plot_box_seq_pair(msim_pl_all,sample_seq,binsz,wr,sample_step,linew,mycc)
ylim([0 0.5])
ylabel('Similarity');box off
set(gca,'position',gcapos);

% rand similarity
subplot(4,2,2);
gcapos = get(gca,'position');
plot_box_seq_pair(msim_rd_all,rand_seq,binsz,wr,sample_step,linew,mycc)
ylim([0 0.5])
box off
set(gca,'position',gcapos);

% plus/minus statistics
subplot(4,2,3);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,1),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
ylabel('Accuracy');box off
set(gca,'position',gcapos);

subplot(4,2,5);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,2),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
ylabel('Precision');box off
set(gca,'position',gcapos);

subplot(4,2,7);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,3),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
ylabel('Recall');box off
xlabel('ensemble%')
set(gca,'position',gcapos);

% rand statistics
subplot(4,2,4);
gcapos = get(gca,'position');
plot_box_seq_single(stats_rd_all(:,:,:,1),rand_seq,wr,sample_step,linew,mycc)
ylim([0 1])
box off
set(gca,'position',gcapos);

subplot(4,2,6);
gcapos = get(gca,'position');
plot_box_seq_single(stats_rd_all(:,:,:,2),rand_seq,wr,sample_step,linew,mycc)
ylim([0 1])
box off
set(gca,'position',gcapos);

subplot(4,2,8);
gcapos = get(gca,'position');
plot_box_seq_single(stats_rd_all(:,:,:,3),rand_seq,wr,sample_step,linew,mycc)
ylim([0 1])
box off
xlabel('total%')
set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee ...
    '_mc_svd_rand_pred_stats_' ge_type '.pdf'])

%% ROC
auc_pl = zeros(length(sample_seq),expt_count);
auc_rd = zeros(length(rand_seq),expt_count);

figure; set(gcf,'color','w','position',[2030 79 974 399])
subplot(1,2,1); hold on
gcapos = get(gca,'position');
step = 1/length(sample_seq);
cc = [interp1(1/64:1/64:1,bkrmap.cmap(:,1),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,2),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,3),step:step:1)];
% cc_light = min(cc+0.3,1);
% cc_dark = max(cc-0.3,0);
cc_dark = cc';
for jj = 1:length(sample_seq)
    for ii = 1:expt_count
        [xx,yy,~,auc_pl(jj,ii)] = perfcurve(repmat(true_label_all{ii},1,num_rand),...
            reshape(squeeze(cell2mat(sim_pl_all(ii,jj,:))),1,[]),1);
%         plot(xx,yy,'color',cc_light(jj,:),'linewidth',linew);
    end
    % average across experiments
    [xx,yy] = perfcurve(repmat(cell2mat(true_label_all),1,num_rand),...
        cell2mat(reshape(sim_pl_all(:,jj,:),1,[])),1);
    if sample_seq(jj)~=0
        plot(xx,yy,'color',cc_dark(jj,:),'linewidth',linew);
    else
        plot(xx,yy,'color',cc_dark(jj,:),'linewidth',3*linew);
    end
end
xlabel('FPR'); ylabel('TPR')
colorbar; colormap(cc_dark); caxis([sample_seq(1) sample_seq(end)])
set(gca,'position',gcapos);

subplot(1,2,2); hold on
gcapos = get(gca,'position');
step = 1/2/length(rand_seq);
cc = [interp1(1/64:1/64:1,bkrmap.cmap(:,1),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,2),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,3),step:step:1)];
% cc_light = min(cc+0.3,1);
% cc_dark = max(cc-0.3,0);
cc_dark = cc';
for jj = 1:length(rand_seq)
    for ii = 1:expt_count
        [xx,yy,~,auc_rd(jj,ii)] = perfcurve(repmat(true_label_all{ii},1,num_rand),...
            reshape(squeeze(cell2mat(sim_rd_all(ii,jj,:))),1,[]),1);
%         plot(xx,yy,'color',cc_light(jj,:),'linewidth',linew);
    end
    % average across experiments
    [xx,yy] = perfcurve(repmat(cell2mat(true_label_all),1,num_rand),...
        cell2mat(reshape(sim_rd_all(:,jj,:),1,[])),1);
    plot(xx,yy,'color',cc_dark(jj,:),'linewidth',linew);
end
xlabel('FPR'); ylabel('TPR')
colorbar; colormap(cc_dark); caxis([rand_seq(1) rand_seq(end)])
set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters','-bestfit',[fig_path 'ensemble_reduction_ROC_curve.pdf'])

%% plot AUC
figure; set(gcf,'color','w','position',[1987 380 659 190])
subplot(1,2,1);hold on;
gcapos = get(gca,'position');
plot_box_seq_single(reshape(auc_pl,1,size(auc_pl,1),size(auc_pl,2)),...
    sample_seq,wr,sample_step,linew,mycc.black)
% for ii = 1:length(xseq)
%     h = boxplot(mean(squeeze(data(:,ii,:)),2),'positions',...
%         xseq(ii),'width',step*wr,'colors',cc);
%     setBoxStyle(h,linew);
% end
% xlim([xseq(1)-step,xseq(end)+step])
% ylim([0 max(data(:))])
% set(gca,'xtick',xseq,'xticklabel',xseq,'XTickLabelRotation',45)
ylim([0 1]); ylabel('AUC')
box off
set(gca,'position',gcapos);

subplot(1,2,2);
gcapos = get(gca,'position');
plot_box_seq_single(reshape(auc_rd,1,size(auc_rd,1),size(auc_rd,2)),...
    rand_seq,wr,sample_step,linew,mycc.black)
ylim([0 1])
box off
xlabel('total%'); ylabel('AUC')
set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters','-bestfit',[fig_path 'ensemble_reduction_AUC.pdf'])

end

