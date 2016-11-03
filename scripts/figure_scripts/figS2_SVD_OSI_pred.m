function [] = figS2_SVD_OSI_pred(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.core;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
OSI_thresh = param.OSI_thresh;
mc_minsz = param.mc_minsz;
linew = param.linew;

load(ccode_path);

%% initialize
num_expt = length(expt_name);
num_crf_osi = zeros(num_expt,1);
num_crf_svd = zeros(num_expt,1);
num_crf = zeros(num_expt,1);
num_osi = zeros(num_expt,1);
num_svd = zeros(num_expt,1);
num_cell = zeros(num_expt,1);

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
stats_pl_all_svd = [];
msim_pl_all_svd = [];
stats_pl_all_osi = [];
msim_pl_all_osi = [];

expt_count = 0;
svd_expt_count = 0;
svd_indx = [];
%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data_high = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1)-num_stim;
    
    %% find ensembles
    % load results: 'core_crf','core_svd'
    load([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat']);
        
    % high OSI
    core_osi = cell(num_stim,1);
    [OSI,OSIstim] = calcOSI(Spikes,vis_stim);
    for ii = 1:num_stim
        core_osi{ii} = find((OSI>OSI_thresh)&(OSIstim==ii));
    end
    
    %% plot ensemble
    rr = 1;
    figure;
    set(gcf,'color','w','position',[2041 290 631 478]);
    set(gcf,'paperpositionmode','auto')
    for ss = 1:num_stim
        
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==ss)';
        
        crf_osi = intersect(core_crf{ss},core_osi{ss});
        num_cell(expt_count) = size(Spikes,1);
        num_crf(expt_count) = length(core_crf{ss});
        num_osi(expt_count) = length(core_osi{ss});
        num_crf_osi(expt_count) = length(crf_osi);
        
        % prediction
        % crf
        [pred.crf{expt_count},cos_sim.crf{expt_count},cos_thresh.crf(expt_count),...
            cos_sim_avg.crf(expt_count,:),acc,prc,rec] = core_cos_sim(core_crf{ss},data_high,true_label);
        pred_stats.crf(expt_count,:) = [acc,prc,rec];
        % osi
        [pred.osi{expt_count},cos_sim.osi{expt_count},cos_thresh.osi(expt_count),...
            cos_sim_avg.osi(expt_count,:),acc,prc,rec] = core_cos_sim(core_osi{ss},data_high,true_label);
        pred_stats.osi(expt_count,:) = [acc,prc,rec];
        
        subplot(2,num_stim,ss);
        plotCoreOverlay(Coord_active,core_crf{ss},core_svd{ss},mycc.orange,...
            mycc.green,rr)
        subplot(2,num_stim,2+ss);
        plotCoreOverlay(Coord_active,core_crf{ss},core_osi{ss},mycc.orange,...
            mycc.gray,rr)
        
        % randomly take out cores - svd
        if ~isempty(core_svd{ss})
            
            svd_expt_count = svd_expt_count+1;
            svd_indx(end+1) = 1;
            
            crf_svd = intersect(core_crf{ss},core_svd{ss});
            num_crf_svd(svd_expt_count) = length(crf_svd);
            num_svd(svd_expt_count) = length(core_svd{ss});
            
            % svd
            [pred.svd{svd_expt_count},cos_sim.svd{svd_expt_count},...
                cos_thresh.svd(svd_expt_count),cos_sim_avg.svd(svd_expt_count,:),...
                acc,prc,rec] = core_cos_sim(core_svd{ss},data_high,true_label);
            pred_stats.svd(svd_expt_count,:) = [acc,prc,rec];
            
            % rand core
            num_core = length(core_svd{ss});
            noncore = setdiff(1:num_node,core_svd{ss});
            core_plus_seq = round(num_core*sample_seq);
            for ii = 1:length(sample_seq)
                for jj = 1:num_rand
                    if core_plus_seq(ii) < 0
                        rand_core = core_osi{ss}(randperm(num_core,num_core+core_plus_seq(ii)));
                    else
                        rand_core = noncore(randperm(length(noncore),core_plus_seq(ii)));
                        rand_core = [core_osi{ss};rand_core'];
                    end
                    % predict
                    [~,~,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data_high,true_label);
                    msim_pl_all_svd(svd_expt_count,ii,jj,:) = sim_avg;
                    stats_pl_all_svd(svd_expt_count,ii,jj,:) = [acc,prc,rec];

                end
            end
        else
            svd_indx(end+1) = 0;
        end
        
        % randomly take out cores - osi
        num_core = length(core_osi{ss});
        noncore = setdiff(1:num_node,core_osi{ss});
        core_plus_seq = round(num_core*sample_seq);
        for ii = 1:length(sample_seq)
            for jj = 1:num_rand
                if core_plus_seq(ii) < 0
                    rand_core = core_osi{ss}(randperm(num_core,num_core+core_plus_seq(ii)));
                else
                    rand_core = noncore(randperm(length(noncore),core_plus_seq(ii)));
                    rand_core = [core_osi{ss};rand_core'];
                end
                
                % predict
                [~,~,~,sim_avg,acc,prc,rec] = core_cos_sim(rand_core,data_high,true_label);
                msim_pl_all_osi(expt_count,ii,jj,:) = sim_avg;
                stats_pl_all_osi(expt_count,ii,jj,:) = [acc,prc,rec];
                
            end
        end
        
    end

    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee '_' ...
        savestr '_svd_osi_core.pdf'])
    
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
h = boxplot(num_crf./num_cell,'positions',stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew);
h = boxplot(num_svd./num_cell(svd_indx==1),'positions',2*stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew);
h = boxplot(num_osi./num_cell,'positions',3*stepsz,'width',ww,'colors',mycc.gray);
setBoxStyle(h,linew);
xlim([0 4*stepsz])
% ylim([0 max(max([num_osi,num_crf,num_svd]./(num_cell*ones(1,3))))])
set(gca,'xcolor','w')
ylabel('cells (%)')
box off

% overlap
subplot(1,2,2); hold on
h = boxplot(num_crf_svd./num_crf(svd_indx==1),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew);
h = boxplot(num_crf_osi./num_crf,'positions',2*stepsz,'width',ww,'colors',mycc.gray);
setBoxStyle(h,linew);
xlim([0 3*stepsz])
ylim([0 1])
ylabel('shared with CRF (%)')
set(gca,'xcolor','w')
box off

saveas(gcf,[fig_path expt_ee '_' savestr '_svd_osi_core_nums.pdf'])

%% plot stats
figure;
set(gcf,'color','w','position',[2038 520 778 221]);
set(gcf,'paperpositionmode','auto')

stepsz = 0.5;
binsz = 0.1;
ww = 0.2;

% mean sim value
subplot(1,4,1); hold on
scatter((stepsz-binsz)*ones(size(cos_sim_avg.crf(:,1))),cos_sim_avg.crf(:,1),...
    30,mycc.orange_light,'+','linewidth',linew);
scatter((stepsz+binsz)*ones(size(cos_sim_avg.crf(:,2))),cos_sim_avg.crf(:,2),...
    30,mycc.orange,'+','linewidth',linew);
scatter((2*stepsz-binsz)*ones(size(cos_sim_avg.svd(:,1))),cos_sim_avg.svd(:,1),...
    30,mycc.green_light,'+','linewidth',linew);
scatter((2*stepsz+binsz)*ones(size(cos_sim_avg.svd(:,2))),cos_sim_avg.svd(:,2),...
    30,mycc.green,'+','linewidth',linew);
scatter((3*stepsz-binsz)*ones(size(cos_sim_avg.osi(:,1))),cos_sim_avg.osi(:,1),...
    30,mycc.gray_light,'+','linewidth',linew);
scatter((3*stepsz+binsz)*ones(size(cos_sim_avg.osi(:,2))),cos_sim_avg.osi(:,2),...
    30,mycc.gray,'+','linewidth',linew);
% crf mean
plot([(stepsz-binsz)*ones(size(cos_sim_avg.crf(:,1))),(stepsz+binsz)*...
    ones(size(cos_sim_avg.crf(:,1)))]',cos_sim_avg.crf','color',mycc.gray);
plot(stepsz-binsz*[1.5 0.5],nanmean(cos_sim_avg.crf(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
plot(stepsz+binsz*[1.5 0.5],nanmean(cos_sim_avg.crf(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
% svd mean
plot([(2*stepsz-binsz)*ones(size(cos_sim_avg.svd(:,1))),(2*stepsz+binsz)*...
    ones(size(cos_sim_avg.svd(:,1)))]',cos_sim_avg.svd','color',mycc.gray);
plot(2*stepsz-binsz*[1.5 0.5],nanmean(cos_sim_avg.svd(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
plot(2*stepsz+binsz*[1.5 0.5],nanmean(cos_sim_avg.svd(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
% osi mean
plot([(3*stepsz-binsz)*ones(size(cos_sim_avg.osi(:,1))),(3*stepsz+binsz)*...
    ones(size(cos_sim_avg.osi(:,1)))]',cos_sim_avg.osi','color',mycc.gray);
plot(3*stepsz-binsz*[1.5 0.5],nanmean(cos_sim_avg.osi(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
plot(3*stepsz+binsz*[1.5 0.5],nanmean(cos_sim_avg.osi(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',3*linew);
xlim([0.2 4*stepsz-0.2])
set(gca,'xcolor','w')
ylabel('Similarity')

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats.crf(:,1),'positions',stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.svd(:,1),'positions',2*stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.osi(:,1),'positions',3*stepsz,'width',ww,'colors',mycc.gray);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
ylabel('accuracy')
set(gca,'xcolor','w')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats.crf(:,2),'positions',stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.svd(:,2),'positions',2*stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.osi(:,2),'positions',3*stepsz,'width',ww,'colors',mycc.gray);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
ylabel('precision')
set(gca,'xcolor','w')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats.crf(:,3),'positions',stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.svd(:,3),'positions',2*stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.osi(:,3),'positions',3*stepsz,'width',ww,'colors',mycc.gray);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
ylabel('recall')
set(gca,'xcolor','w')
box off

saveas(gcf,[fig_path expt_ee '_svd_osi_core_pred_' ge_type '_stats.pdf'])

%% plot reduction results
binsz = 0.02;
wr = 0.3;

figure
set(gcf,'color','w')
set(gcf,'position',[1991 632 1016 298]);
set(gcf,'paperpositionmode','auto')

% ------------------ svd ------------------- %
% plus/minus similarity
subplot(2,4,1);
gcapos = get(gca,'position');
plot_box_seq_pair(msim_pl_all_svd,sample_seq,binsz,wr,sample_step,linew,...
    mycc.green_light,mycc.green)
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylim([0 0.5])
ylabel('similarity');box off
set(gca,'position',gcapos);

% plus/minus statistics
subplot(2,4,2);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_svd(:,:,:,1),sample_seq,wr,sample_step,linew,mycc.green)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('accuracy');box off
set(gca,'position',gcapos);

subplot(2,4,3);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_svd(:,:,:,2),sample_seq,wr,sample_step,linew,mycc.green)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('precision');box off
set(gca,'position',gcapos);

subplot(2,4,4);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_svd(:,:,:,3),sample_seq,wr,sample_step,linew,mycc.green)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('recall');box off
set(gca,'position',gcapos);

% ------------------ osi ------------------- %
% plus/minus similarity
subplot(2,4,5);
gcapos = get(gca,'position');
plot_box_seq_pair(msim_pl_all_osi,sample_seq,binsz,wr,sample_step,linew,...
    mycc.gray_light,mycc.gray)
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylim([0 0.5])
ylabel('similarity');box off
xlabel('core neurons (%)')
set(gca,'position',gcapos);

% plus/minus statistics
subplot(2,4,6);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_osi(:,:,:,1),sample_seq,wr,sample_step,linew,mycc.gray)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('accuracy');box off
xlabel('core neurons (%)')
set(gca,'position',gcapos);

subplot(2,4,7);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_osi(:,:,:,2),sample_seq,wr,sample_step,linew,mycc.gray)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('precision');box off
set(gca,'position',gcapos);
xlabel('core neurons (%)')

subplot(2,4,8);
gcapos = get(gca,'position');
plot_box_seq_single(stats_pl_all_osi(:,:,:,3),sample_seq,wr,sample_step,linew,mycc.gray)
ylim([0 1])
set(gca,'xtick',sample_seq,'xticklabel',num2str(100*(1+sample_seq')));
ylabel('recall');box off
xlabel('core neurons (%)')
set(gca,'position',gcapos);

print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_ee ...
    '_SVD_OSI_reduction_pred_stats_' ge_type '.pdf'])

end