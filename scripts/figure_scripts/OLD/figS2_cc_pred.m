function [] = figS2_cc_pred(param)

% parse parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.cc;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
savestr = param.savestr;
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

sample_step = 0.1;
num_rand = 100;

sample_seq = -0.9:sample_step:0;
[~,indx] = min(abs(sample_seq));
sample_seq(indx) = 0;

% initialize
stats_pl_all = [];
msim_pl_all = [];

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
    core_crf = find_core_max_clique(best_model.graph,num_stim,mc_minsz);
    
    %% plot ep sum
%     figure; set(gcf,'color','w','position',[2154 560 479 377])
%     cc_range = [];
%     ep_range = [-1.5 0.1];
%     plotGraphModelHighlightEP(cc_graph,coords,model.ep_on,...
%         cc_range,ep_range,gmap.cmap,[]);
    
    %% plot ensemble
    rr = 1;
    hf = figure;
    set(gcf,'color','w','position',[2041 430 543 338]);
    set(gcf,'paperpositionmode','auto')
    for ss = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ss)';
        
        crf_cc = intersect(core_crf{ss},core_cc{ss});
        num_cell(expt_count) = size(Spikes,1);
        num_crf(expt_count) = length(core_crf{ss});
        num_cc(expt_count) = length(core_cc{ss});
        num_crf_cc(expt_count) = length(crf_cc);
        
        % prediction
        [pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
            cos_sim_avg(expt_count,:),acc,prc,rec] = core_cos_sim(core_cc{ss},data_high,true_label);
        pred_stats(expt_count,:) = [acc,prc,rec];
        
        subplot(1,num_stim,ss);
        plotCoreOverlay(Coord_active,core_crf{ss},core_cc{ss},mycc.orange,...
            mycc.gray,rr)
        
        % randomly take out cores
        num_core = length(core_cc{ss});
        noncore = setdiff(1:num_node,core_cc{ss});
        core_plus_seq = round(num_core*sample_seq);
        for ii = 1:length(sample_seq)
            for jj = 1:num_rand
                if core_plus_seq(ii) < 0
                    rand_core = core_cc{ss}(randperm(num_core,num_core+core_plus_seq(ii)));
                else
                    rand_core = noncore(randperm(length(noncore),core_plus_seq(ii)));
                    rand_core = [core_cc{ss};rand_core'];
                end
                
                % predict
                [~,~,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data_high,true_label);
                msim_pl_all(expt_count,ii,jj,:) = sim_avg;
                stats_pl_all(expt_count,ii,jj,:) = [acc,prc,rec];
                
            end
        end
        
        % cosine similarity
        true_label = double(vis_stim_high==ss)';

        [pred{expt_count},cos_sim{expt_count},cos_thresh(expt_count),...
            cos_sim_avg(expt_count,:),acc,prc,rec] = core_cos_sim(core_cc{ss},data_high,true_label);
        pred_stats(expt_count,:) = [acc,prc,rec];
        
    end

    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee '_' ...
        savestr '_mc_cc_core.pdf'])
    
end

%% plot core numbers
stepsz = 0.5;
ww = 0.3;

figure
set(gcf,'color','w');
set(gcf,'position',[2096 452 447 233])
set(gcf,'paperpositionmode','auto')

% ensemble number
subplot(1,2,1); hold on
h = boxplot(num_cc./num_cell,'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew);
xlim([0 2*stepsz])
ylim([0 max(num_cc./num_cell)])
set(gca,'xcolor','w')
ylabel('cells (%)')
box off

% CRF+osi
subplot(1,2,2); hold on
h = boxplot(num_crf_cc./num_crf,'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew);
xlim([0 3*stepsz])
ylim([0 1])
ylabel('CC in CRF (%)')
set(gca,'xcolor','w')
box off

saveas(gcf,[fig_path expt_ee '_' savestr '_mc_cc_core_nums.pdf'])

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

%% reduction
binsz = 0.02;
wr = 0.3;

figure
set(gcf,'color','w')
set(gcf,'position',[1991 327 650 333]);
set(gcf,'paperpositionmode','auto')

% plus/minus similarity
subplot(2,2,1);
gcapos = get(gca,'position');
plot_box_seq_pair(msim_pl_all,sample_seq,binsz,wr,sample_step,linew,mycc)
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylim([0 0.6])
ylabel('Similarity');box off
set(gca,'position',gcapos);

% plus/minus statistics
subplot(2,2,2);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,1),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('Accuracy');box off
set(gca,'position',gcapos);

subplot(2,2,3);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,2),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('Precision');box off
set(gca,'position',gcapos);
xlabel('core neurons (%)')

subplot(2,2,4);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all(:,:,:,3),sample_seq,wr,sample_step,linew,mycc)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('Recall');box off
xlabel('core neurons (%)')
set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee ...
    '_CC_reduction_pred_stats_' ge_type '.pdf'])

end