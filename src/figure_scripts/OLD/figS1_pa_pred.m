function [] = figS1_pa_pred(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.pa;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
OSI_thresh = param.OSI_thresh;
mc_minsz = param.mc_minsz;
linew = param.linew;

load(ccode_path);

%% initialize
num_expt = length(expt_name);
num_crf = zeros(num_expt,1);
num_cell = zeros(num_expt,1);

cos_sim = {};
cos_sim_avg = [];
cos_thresh = [];
pred = {};
pred_stats = [];

expt_count = 0;

%% graph properties
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data_high = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    num_stim = length(setdiff(unique(vis_stim),0));
    
    %% find ensembles
    core_crf = find_core_max_clique(best_model.graph,num_stim,mc_minsz);
    
    %% plot ensemble
    rr = 1;
    figure;
    set(gcf,'color','w','position',[2041 430 543 338]);
    set(gcf,'paperpositionmode','auto')
    for ii = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ii)';
        
        num_crf(expt_count) = length(core_crf{ii});
        num_cell(expt_count) = size(best_model.graph,1)-num_stim;
        
        % prediction
        [pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
            cos_sim_avg(expt_count,:),acc,prc,rec] = core_cos_sim(core_crf{ii},...
            data_high,true_label');
        pred_stats(expt_count,:) = [acc,prc,rec];
        
%         subplot(1,num_stim,ii);
%         plotCoreOverlay(Coord_active,core_crf{ii},core_osi{ii},mycc.orange,...
%             mycc.gray,rr)
        
    end

%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee '_' ...
%         savestr '_mc_osi_core.pdf'])
    
    
end

%% plot stats
stepsz = 0.5;
ww = 0.3;

figure
set(gcf,'color','w');
set(gcf,'position',[2096 452 447 233])
set(gcf,'paperpositionmode','auto')
hold on
h = boxplot(num_crf./num_cell,'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew);
xlim([0 2*stepsz])
ylim([0 max(num_crf./num_cell)])
set(gca,'xcolor','w')
ylabel('cells (%)')
box off

saveas(gcf,[fig_path 'pa_' savestr '_crf_core_nums.pdf'])

%% plot stats
figure;
set(gcf,'color','w','position',[2038 520 778 221]);
set(gcf,'paperpositionmode','auto')

stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

% mean sim value
subplot(1,4,1); hold on
scatter((stepsz-binsz)*ones(size(cos_sim_avg(:,1))),cos_sim_avg(:,1),30,mycc.gray,'+','linewidth',linew);
scatter((stepsz+binsz)*ones(size(cos_sim_avg(:,2))),cos_sim_avg(:,2),30,mycc.black,'+','linewidth',linew);
plot([(stepsz-binsz)*ones(size(cos_sim_avg(:,1))),(stepsz+binsz)*ones(size(cos_sim_avg(:,1)))]',...
    cos_sim_avg','color',mycc.gray);
plot([stepsz-binsz*1.5,stepsz-binsz*0.5],nanmean(cos_sim_avg(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
plot([stepsz+binsz*0.5,stepsz+binsz*1.5],nanmean(cos_sim_avg(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
xlim([0.2 2*stepsz-0.2])
set(gca,'xcolor','w')
ylabel('Similarity')

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats(:,1),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
ylabel('Accuracy')
set(gca,'xcolor','w')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats(:,2),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
ylabel('Precision')
set(gca,'xcolor','w')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats(:,3),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
ylabel('Recall')
set(gca,'xcolor','w')
box off

saveas(gcf,[fig_path 'pa_crf_core_pred_' ge_type '_stats.pdf'])

end
