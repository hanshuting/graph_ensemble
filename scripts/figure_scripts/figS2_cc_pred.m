function [] = figS2_cc_pred(param)

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
mc_minsz = param.mc_minsz;

load(ccode_path);
load(rwbmap);

%%
cos_sim = {};
cos_sim_avg = [];
cos_thresh = [];
pred = {};
pred_stats = [];

expt_count = 0;
for n = 1:length(expt_name)
    expt_ee = ee{n}{1};
    load([param.data_path expt_name{n} '\' expt_name{n} '.mat']);
    fprintf('Processing %s_%s...\n',expt_name{n},expt_ee);
    
    % set path
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    cc_path = [result_path_base '\' expt_name{n} '\cc\']; 
    
    % load data
    best_model = load([model_path  expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    load([cc_path  expt_name{n} '_' expt_ee '_cc_graph.mat']);    
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data_high = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    num_stim = length(unique(vis_stim))-1;
    
    crf_graph = best_model.graph;
    edge_pot = best_model.edge_pot;
    num_node = size(crf_graph,1);

    core_cc = find_core_max_clique(cc_graph,num_stim,mc_minsz);
        
    %% cosine similarity
    for ii = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ii)';

        [pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
            cos_sim_avg(expt_count,:),acc,prc,rec] = core_cos_sim(core_cc{ii},data_high,true_label);
        pred_stats(expt_count,:) = [acc,prc,rec];
        
    end
    
end

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
ylabel('Similarity')

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats(:,1),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
ylabel('Accuracy')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats(:,2),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
ylabel('Precision')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats(:,3),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew)
xlim([0 2*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2]*stepsz);
ylabel('Recall')
box off

saveas(gcf,[fig_path expt_ee '_mc_cc_core_pred_' ge_type '_stats.pdf'])

end