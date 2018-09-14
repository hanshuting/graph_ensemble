function [] = fig5_plot_crf_pred_cos_CRFSVD(param)

% parse parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.stats;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
linew = param.linew;

load(ccode_path);
load(rwbmap);

%% initialize
cos_sim = struct();
cos_sim_avg = struct();
cos_thresh = struct();
pred = struct();
pred_stats = struct();

%% graph properties
expt_count = 0;
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_node = size(best_model.graph,1);
    num_stim = length(unique(vis_stim))-1;
    
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    load([result_path_base '\' expt_name{n} '\core\' expt_ee '_mc_svd_core_' ...
        ge_type '.mat']);
    
    %% cosine similarity
    for ii = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ii)';
        
        % SVD
        [pred.svd{expt_count},cos_sim.svd{expt_count},cos_thresh.svd{expt_count},...
            cos_sim_avg.svd{expt_count},acc,prc,rec] = core_cos_sim(core_svd{ii},data,true_label);
        pred_stats.svd(expt_count,:) = [acc,prc,rec];
        
        % CRF
        [pred.crf{expt_count},cos_sim.crf{expt_count},cos_thresh.crf{expt_count},...
            cos_sim_avg.crf{expt_count},acc,prc,rec] = core_cos_sim(core_crf{ii},data,true_label);
        pred_stats.crf(expt_count,:) = [acc,prc,rec];
        
    end
    
    %% plot cos sim
    pind = expt_count-num_stim+1:expt_count;
    sind = [1,3,2,4]; % sort by vis stim type
    pred_mat = reshape([cell2mat(pred.crf(pind));cell2mat(pred.svd(pind))]',[],num_stim*2)';
    plot_pred_raster(pred_mat(sind,:),vis_stim_high,cmap)
    
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
        expt_ee '_mc_svd_core_pred_raster_' ge_type '.pdf'])
    
end

%% plot stats
figure;
set(gcf,'color','w','position',[2041 533 993 235]);
set(gcf,'paperpositionmode','auto')

stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

% mean sim value
subplot(1,4,1); hold on
% svd model
mcs = cell2mat(cos_sim_avg.svd');
scatter((stepsz-binsz)*ones(size(mcs(:,1))),mcs(:,1),30,mycc.green_light,'+','linewidth',linew);
scatter((stepsz+binsz)*ones(size(mcs(:,2))),mcs(:,2),30,mycc.green,'+','linewidth',linew);
plot([(stepsz-binsz)*ones(size(mcs(:,1))),(stepsz+binsz)*ones(size(mcs(:,1)))]',...
    mcs','color',mycc.gray);
plot([stepsz-binsz*1.5,stepsz-binsz*0.5],nanmean(mcs(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
plot([stepsz+binsz*0.5,stepsz+binsz*1.5],nanmean(mcs(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
% crf model
mcs = cell2mat(cos_sim_avg.crf');
scatter((2*stepsz-binsz)*ones(size(mcs(:,1))),mcs(:,1),30,mycc.orange_light,'+','linewidth',linew);
scatter((2*stepsz+binsz)*ones(size(mcs(:,2))),mcs(:,2),30,mycc.orange,'+','linewidth',linew);
plot([(2*stepsz-binsz)*ones(size(mcs(:,1))),(2*stepsz+binsz)*ones(size(mcs(:,1)))]',...
    mcs','color',mycc.gray);
plot([2*stepsz-binsz*1.5,2*stepsz-binsz*0.5],mean(mcs(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
plot([2*stepsz+binsz*0.5,2*stepsz+binsz*1.5],mean(mcs(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
xlim([0.2 3*stepsz-0.2])
% ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
set(gca,'xticklabel',{'SVD','CRF'})
ylabel('Similarity')

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats.svd(:,1),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,1),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
set(gca,'xticklabel',{'SVD','CRF'})
ylabel('Accuracy')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats.svd(:,2),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,2),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
set(gca,'xticklabel',{'SVD','CRF'})
ylabel('Precision')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats.svd(:,3),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,3),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
xlim([0 3*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
set(gca,'xticklabel',{'SVD','CRF'})
ylabel('Recall')
box off

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee '_mc_svd_core_pred_'...
     ge_type '_stats.pdf'])

end
