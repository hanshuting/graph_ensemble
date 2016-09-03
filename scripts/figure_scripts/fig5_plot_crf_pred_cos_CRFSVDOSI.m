function [] = fig5_plot_crf_pred_cos_CRFSVDOSI(param)

% parse parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
k = param.k;
p = param.p;
cc_type = param.cc_type;
comm_type = param.comm_type;
data_path = param.data_path;
fig_path = param.fig_path.stats;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;

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

    model_path = [result_path_base expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee '_loopy_best_model.mat']);
    num_node = size(best_model.graph,1);
    num_stim = length(unique(vis_stim))-1;
    
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    load([result_path_base expt_name{n} '\core\' expt_ee '_mc_core.mat']);
    
    %% cosine similarity
    for ii = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ii)';
        
        % SVD
        [pred.svd{expt_count},cos_sim.svd{expt_count},cos_thresh.svd{expt_count},...
            cos_sim_avg.svd{expt_count},acc,prc,rec] = core_cos_sim(core_svd{ii},data,true_label);
        pred_stats.svd(expt_count,:) = [acc,prc,rec];
        
        % OSI
        [pred.osi{expt_count},cos_sim.osi{expt_count},cos_thresh.osi{expt_count},...
            cos_sim_avg.osi{expt_count},acc,prc,rec] = core_cos_sim(core_osi{ii},data,true_label);
        pred_stats.osi(expt_count,:) = [acc,prc,rec];
        
        % CRF
        [pred.crf{expt_count},cos_sim.crf{expt_count},cos_thresh.crf{expt_count},...
            cos_sim_avg.crf{expt_count},acc,prc,rec] = core_cos_sim(core_crf{ii},data,true_label);
        pred_stats.crf(expt_count,:) = [acc,prc,rec];
        
    end
    
    %% plot cos sim
    pind = expt_count-num_stim+1:expt_count;
    sind = [1,3,5,2,4,6]; % sort by vis stim type
    pred_mat = reshape([cell2mat(pred.crf(pind));cell2mat(pred.svd(pind));...
        cell2mat(pred.osi(pind))]',[],num_stim*3)';
    sim_mat = reshape([cell2mat(cos_sim.crf(pind));cell2mat(cos_sim.svd(pind));...
        cell2mat(cos_sim.osi(pind))]',[],num_stim*3)';
    thresh_mat = [cell2mat(cos_thresh.crf(pind)),cell2mat(cos_thresh.svd(pind)),...
        cell2mat(cos_thresh.osi(pind))]';
    plot_pred(pred_mat(sind,:),sim_mat(sind,:),thresh_mat(sind),vis_stim_high,cmap)
    
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' expt_ee '_mc_core_pred.pdf'])
    
end

%% plot stats
figure;
set(gcf,'color','w','position',[2041 533 1195 235]);
set(gcf,'paperpositionmode','auto')

stepsz = 1;
binsz = 0.1;
ww = 0.3;

% mean sim value
subplot(1,4,1); hold on
mcs = cell2mat(cos_sim_avg.crf');
h = boxplot(mcs(:,1),'positions',stepsz-binsz,'width',ww,'colors',mycc.red_light);
set(h(7,:),'visible','off')
h = boxplot(mcs(:,2),'positions',stepsz+binsz,'width',ww,'colors',mycc.red);
set(h(7,:),'visible','off')
mcs = cell2mat(cos_sim_avg.svd');
h = boxplot(mcs(:,1),'positions',2*stepsz-binsz,'width',ww,'colors',mycc.blue_light);
set(h(7,:),'visible','off')
h = boxplot(mcs(:,2),'positions',2*stepsz+binsz,'width',ww,'colors',mycc.blue);
set(h(7,:),'visible','off')
mcs = cell2mat(cos_sim_avg.osi');
h = boxplot(mcs(:,1),'positions',3*stepsz-binsz,'width',ww,'colors',mycc.purple_light);
set(h(7,:),'visible','off')
h = boxplot(mcs(:,2),'positions',3*stepsz+binsz,'width',ww,'colors',mycc.purple);
set(h(7,:),'visible','off')
xlim([0 4*stepsz])
ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'CRF','SVD','OSI'})
ylabel('Similarity')

% accuracy
stepsz = 0.5;
ww = 0.2;
subplot(1,4,2); hold on
h = boxplot(pred_stats.crf(:,1),'positions',stepsz,'width',ww,'colors',mycc.red);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.svd(:,1),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.osi(:,1),'positions',3*stepsz,'width',ww,'colors',mycc.purple);
set(h(7,:),'visible','off')
xlim([0 4*stepsz])
ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'CRF','SVD','OSI'})
ylabel('Accuracy')

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats.crf(:,2),'positions',stepsz,'width',ww,'colors',mycc.red);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.svd(:,2),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.osi(:,2),'positions',3*stepsz,'width',ww,'colors',mycc.purple);
set(h(7,:),'visible','off')
xlim([0 4*stepsz])
ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'CRF','SVD','OSI'})
ylabel('Precision')

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats.crf(:,3),'positions',stepsz,'width',ww,'colors',mycc.red);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.svd(:,3),'positions',2*stepsz,'width',ww,'colors',mycc.blue);
set(h(7,:),'visible','off')
h = boxplot(pred_stats.osi(:,3),'positions',3*stepsz,'width',ww,'colors',mycc.purple);
set(h(7,:),'visible','off')
xlim([0 4*stepsz])
ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'CRF','SVD','OSI'})
ylabel('Recall')

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee '_mc_core_pred_stats.pdf'])

end
