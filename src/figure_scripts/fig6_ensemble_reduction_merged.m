function [] = fig6_ensemble_reduction_merged(param)

% parse parameters
expt_name = param.expt_name;
ee = param.ee;
data_path = param.data_path;
fig_path = param.fig_path.stats;
ge_type = param.ge_type;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
linew = param.linew;
p = param.p;

bkrmap = load(param.bkrmap);
load(ccode_path);

sample_step = 0.1;
num_rand = 100;

%% ensemble predictions
num_expt = length(expt_name);
sample_seq = -0.9:sample_step:0.9;
[~,indx] = min(abs(sample_seq));
sample_seq(indx) = 0;
% rand_seq = 0.1:sample_step:1;
rand_seq = sample_seq;

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
    num_stim = length(unique(setdiff(vis_stim,0)));
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
                num_rd_core = length(rand_core);
                
                % predict
                [~,sim_core,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data,true_label);
                sim_pl_all{expt_count,ii,jj} = sim_core;
                msim_pl_all(expt_count,ii,jj,:) = sim_avg;
                stats_pl_all(expt_count,ii,jj,:) = [acc,prc,rec];
                
                % randomly sample neurons
                rand_core = randperm(num_node,num_rd_core);
                % predict
                [~,sim_core,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data,true_label);
                sim_rd_all{expt_count,ii,jj} = sim_core;
                msim_rd_all(expt_count,ii,jj,:) = sim_avg;
                stats_rd_all(expt_count,ii,jj,:) = [acc,prc,rec];
                
            end
        end
    end
    
end

%% significance test
pval_pl = struct();
pval_rd = struct();
full_indx = 1/sample_step;

% similarity
for ii = 1:length(sample_seq)
    pval.sim(ii) = ranksum(mean(squeeze(msim_pl_all(:,ii,:,2)),2),...
        mean(squeeze(msim_pl_all(:,full_indx,:,2)),2));
end
for ii = 1:length(rand_seq)
    pval_rd.sim(ii) = ranksum(mean(squeeze(msim_rd_all(:,ii,:,2)),2),...
        mean(squeeze(msim_pl_all(:,full_indx,:,2)),2));
end

% accuracy
for ii = 1:length(sample_seq)
    pval_pl.acc(ii) = ranksum(mean(squeeze(stats_pl_all(:,ii,:,1)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,1)),2));
end
for ii = 1:length(rand_seq)
    pval_rd.acc(ii) = ranksum(mean(squeeze(stats_rd_all(:,ii,:,1)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,1)),2));
end

% precision
for ii = 1:length(sample_seq)
    pval_pl.prc(ii) = ranksum(mean(squeeze(stats_pl_all(:,ii,:,2)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,2)),2));
end
for ii = 1:length(rand_seq)
    pval_rd.prc(ii) = ranksum(mean(squeeze(stats_rd_all(:,ii,:,2)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,2)),2));
end

% recall
for ii = 1:length(sample_seq)
    pval_pl.rec(ii) = ranksum(mean(squeeze(stats_pl_all(:,ii,:,3)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,3)),2));
end
for ii = 1:length(rand_seq)
    pval_rd.rec(ii) = ranksum(mean(squeeze(stats_rd_all(:,ii,:,3)),2),...
        mean(squeeze(stats_pl_all(:,full_indx,:,3)),2));
end

%% ROC
auc_pl = zeros(length(sample_seq),expt_count);
auc_rd = zeros(length(rand_seq),expt_count);

xvec = 0:0.02:1;

figure; set(gcf,'position',[2028 485 357 253],'color','w');hold on
% gcapos = get(gca,'position');
step = 1/length(sample_seq);

cc_pl = [interp1(1/64:1/64:1,bkrmap.cmap(:,1),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,2),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,3),step:step:1)]';

% cc_rd = [interp1(1/64:1/64:1,bkrmap.cmap(:,1),step:step:1);...
%     interp1(1/64:1/64:1,bkrmap.cmap(:,2),step:step:1);...
%     interp1(1/64:1/64:1,bkrmap.cmap(:,3),step:step:1)]';
cc_rd = [interp1(1/64:1/64:1,bkrmap.cmap(:,1),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,2),step:step:1);...
    interp1(1/64:1/64:1,bkrmap.cmap(:,3),step:step:1)]';

xx_pl = cell(expt_count,1);
yy_pl = cell(expt_count,1);
xx_rd = cell(expt_count,1);
yy_rd = cell(expt_count,1);
for jj = 1:length(sample_seq)
    
    for ii = 1:expt_count
        [xx_pl{ii},yy_pl{ii},~,auc_pl(jj,ii)] = perfcurve(repmat(reshape(true_label_all{ii},1,[]),1,num_rand),...
            reshape(squeeze(cell2mat(sim_pl_all(ii,jj,:))),1,[]),1);
        [xx_rd{ii},yy_rd{ii},~,auc_rd(jj,ii)] = perfcurve(repmat(reshape(true_label_all{ii},1,[]),1,num_rand),...
            reshape(squeeze(cell2mat(sim_rd_all(ii,jj,:))),1,[]),1);
    end
    
    % average random
    ymat = zeros(expt_count,length(xvec));
    for ii = 1:expt_count
        [~,uid] = unique(xx_rd{ii,1});
        ymat(ii,:) = interp1(xx_rd{ii}(uid),yy_rd{ii}(uid),xvec);
    end
    if sample_seq(jj)~=0
        plot(xvec,mean(ymat,1),'color',cc_rd(jj,:),'linewidth',linew);
    else
        plot(xvec,mean(ymat,1),'color',cc_rd(jj,:),'linewidth',3*linew);
    end
    
    % average plus/minus
    ymat = zeros(expt_count,length(xvec));
    for ii = 1:expt_count
        [~,uid] = unique(xx_pl{ii,1});
        ymat(ii,:) = interp1(xx_pl{ii}(uid),yy_pl{ii}(uid),xvec);
    end
    if sample_seq(jj)~=0
        plot(xvec,mean(ymat,1),'color',cc_pl(jj,:),'linewidth',linew);
    else
        plot(xvec,mean(ymat,1),'color',cc_pl(jj,:),'linewidth',3*linew);
    end
    
end
xlabel('FPR'); ylabel('TPR')
colorbar; colormap(cc_rd); caxis([sample_seq(1) sample_seq(end)])
% set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters',[fig_path expt_ee ...
    '_core_rand_pred_ROC_' ge_type '_merged.pdf'])


%% plot 
binsz = 0.02;
wr = 0.3;

hf = figure;
set(gcf,'color','w')
set(gcf,'position',[2041 361 667 612]);
set(gcf,'paperpositionmode','auto')

% similarity
subplot(3,2,1);
gcapos = get(gca,'position');
plot_box_seq_pair(cat(4,msim_rd_all(:,:,:,2),msim_pl_all(:,:,:,2)),...
    sample_seq,binsz,wr,sample_step,linew,mycc.black,mycc.orange)
ylim([0 0.5])
ylabel('Similarity');box off
set(gca,'position',gcapos);
for ii = 1:length(sample_seq)
    if pval_pl.sim(ii)<p
        scatter(sample_seq(ii),0.5,'k*');
    end
end

subplot(3,2,2);
gcapos = get(gca,'position');
plot_box_seq_pair(cat(4,stats_rd_all(:,:,:,1),stats_pl_all(:,:,:,1)),...
    sample_seq,binsz,wr,sample_step,linew,mycc.black,mycc.orange)
ylim([0 1])
ylabel('Accuracy');box off
set(gca,'position',gcapos);
for ii = 1:length(sample_seq)
    if pval_pl.acc(ii)<p
        scatter(sample_seq(ii),1,'k*');
    end
end

subplot(3,2,4);
gcapos = get(gca,'position');
plot_box_seq_pair(cat(4,stats_rd_all(:,:,:,2),stats_pl_all(:,:,:,2)),...
    sample_seq,binsz,wr,sample_step,linew,mycc.black,mycc.orange)
ylim([0 1])
ylabel('Precision');box off
set(gca,'position',gcapos);
for ii = 1:length(sample_seq)
    if pval_pl.prc(ii)<p
        scatter(sample_seq(ii),1,'k*');
    end
end

subplot(3,2,6);
gcapos = get(gca,'position');
plot_box_seq_pair(cat(4,stats_rd_all(:,:,:,3),stats_pl_all(:,:,:,3)),...
    sample_seq,binsz,wr,sample_step,linew,mycc.black,mycc.orange)
ylim([0 1])
ylabel('Recall');box off
set(gca,'position',gcapos);
for ii = 1:length(sample_seq)
    if pval_pl.rec(ii)<p
        scatter(sample_seq(ii),1,'k*');
    end
end

% --------- plot AUC -----------
% significance test
for ii = 1:length(sample_seq)
    pval_pl.auc(ii) = ranksum(auc_pl(ii,:),auc_pl(full_indx,:));
end
for ii = 1:length(rand_seq)
    pval_rd.auc(ii) = ranksum(auc_rd(ii,:),auc_pl(full_indx,:));
end

% plot
subplot(3,2,3);hold on;
gcapos = get(gca,'position'); % ############# FIX HERE
plot_box_seq_single(reshape(auc_rd,1,size(auc_rd,1),size(auc_rd,2)),...
    sample_seq-0.05,wr,sample_step,linew,mycc.black)
hold on
plot_box_seq_single(reshape(auc_pl,1,size(auc_pl,1),size(auc_pl,2)),...
    sample_seq+0.05,wr,sample_step,linew,mycc.orange)
ylim([0 1]); ylabel('AUC')
box off; xlabel('core %')
set(gca,'position',gcapos);
for ii = 1:length(sample_seq)
    if pval_pl.auc(ii)<p
        scatter(sample_seq(ii),1,'k*');
    end
end

print(hf,'-dpdf','-painters','-bestfit',[fig_path expt_ee ...
    '_core_rand_pred_stats_' ge_type '.pdf'])

end

