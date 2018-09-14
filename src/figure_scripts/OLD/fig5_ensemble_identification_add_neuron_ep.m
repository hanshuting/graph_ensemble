function [] = fig5_ensemble_identification_add_neuron_ep(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.core;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
OSI_thresh = param.OSI_thresh;
mc_minsz = param.mc_minsz;
linew = param.linew;

load(ccode_path);

%% initialize
num_expt = length(expt_name);
num_crf_svd = zeros(num_expt,1);
num_crf_osi = zeros(num_expt,1);
num_crf = zeros(num_expt,1);
num_svd = zeros(num_expt,1);
num_osi = zeros(num_expt,1);
num_cell = zeros(num_expt,1);

cos_sim = struct();
cos_sim_avg = struct();
cos_thresh = struct();
pred = struct();
pred_stats = struct();

num_mc_in_core = [];

expt_count = 0;

%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G)';
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1)-num_stim;
    
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    %% find ensembles
    % this code now only works with two stimuli
    % SVD
    core_svd = cell(num_stim,1);
    for ii = 1:num_stim
        for jj = 1:length(svd_data.svd_state)
            if strcmp(num2str(ii),svd_data.svd_state{jj})
                core_svd{ii} = svd_data.core_svd{jj};
                break;
            end
        end
    end
    
    % edge potential
    perc_ep = 0.3;
    core_ep = cell(num_stim,1);
    for ii = 1:num_stim
        ens = setdiff(find(best_model.graph(num_node+ii,:)~=0),num_node+ii);
        ep_sum = sum(best_model.edge_pot(ens,ens),2);
        core_ep{ii} = ens(ep_sum>quantile(ep_sum,1-perc_ep));
    end
    
%     figure;set(gcf,'color','w')
%     for ii = 1:num_stim
%         ens = setdiff(find(best_model.graph(num_node+ii,:)~=0),num_node+ii);
%         subplot(1,num_stim,ii)
%         plotGraphModelHighlightEP(best_model.graph(ens,ens),...
%             Coord_active(ens,:),best_model.ep_on(ens,ens),[]);
%     end
%     suptitle(expt_name{n})
    
    % CRF max clique
    [core_crf,mc_in_core] = find_core_max_clique(best_model.graph,num_stim,mc_minsz);
    num_in_core = cellfun('length',mc_in_core);
    num_mc_in_core(end+1:end+length(num_in_core)) = num_in_core;
    
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
        
        % ep
        [pred.ep{expt_count},cos_sim.ep{expt_count},cos_thresh.ep{expt_count},...
            cos_sim_avg.ep{expt_count},acc,prc,rec] = core_cos_sim(core_ep{ii},data,true_label);
        pred_stats.ep(expt_count,:) = [acc,prc,rec];
        
    end
    
    %% plot ensemble
%     rr = 1;
%     figure;
%     set(gcf,'color','w','position',[2041 430 543 338]);
%     set(gcf,'paperpositionmode','auto')
%     for ii = 1:num_stim
%         
%         expt_count = expt_count+1;
%         
%         crf_svd = intersect(core_crf{ii},core_svd{ii});
%         
%         num_cell(expt_count) = size(Spikes,1);
%         num_crf(expt_count) = length(core_crf{ii});
%         num_svd(expt_count) = length(core_svd{ii});
%         num_crf_svd(expt_count) = length(crf_svd);
%         
%         subplot(1,num_stim,ii);
%         plotCoreOverlay(Coord_active,core_crf{ii},core_svd{ii},mycc.orange,...
%             mycc.green,rr)
%         
%     end
%     
%     print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee '_mc_svd_core.pdf'])
%     
    %% save result
%     save([result_path_base '\' expt_name{n} '\core\' expt_ee '_mc_svd_core_' ...
%         ge_type '.mat'],'core_svd','core_crf','-v7.3');
    
end

%% plot number of maximal cliques in core
% figure;set(gcf,'color','w','position',[2008 525 577 219])
% boxwd = 0.2;hold on;
% h = boxplot(num_mc_in_core,'positions',0.5,'width',...
%     boxwd,'colors',mycc.black);
% setBoxStyle(h,linew);
% xlim([0 1])
% gcapos = get(gca,'position');
% ylabel('# maximal cliques in core')
% set(gca,'linewidth',linew)
% set(gca,'position',gcapos);
% box off
% 
% saveas(gcf,[fig_path 'num_mc_in_core.pdf']);

%% plot stats
% stepsz = 0.5;
% ww = 0.3;
% 
% figure
% set(gcf,'color','w');
% set(gcf,'position',[2096 452 447 233])
% set(gcf,'paperpositionmode','auto')
% 
% % ensemble number
% subplot(1,2,1); hold on
% h = boxplot(num_svd./num_cell,'positions',stepsz,'width',ww,'colors',mycc.green);
% setBoxStyle(h,linew);
% h = boxplot(num_crf./num_cell,'positions',2*stepsz,'width',ww,'colors',mycc.orange);
% setBoxStyle(h,linew);
% xlim([0 3*stepsz])
% ylim([0 max([num_svd./num_cell;num_crf./num_cell])])
% set(gca,'xtick',[1,2]*stepsz,'xticklabel',{'SVD','CRF'})
% ylabel('Cells %')
% box off
% 
% % CRF+SVD
% subplot(1,2,2); hold on
% h = boxplot(num_crf_svd./(num_crf+num_svd),'positions',stepsz,'width',ww,'colors',mycc.black);
% setBoxStyle(h,linew);
% xlim([0 3*stepsz])
% ylim([0 1])
% ylabel('Nshared%')
% set(gca,'xcolor','w')
% box off
% 
% print(gcf,'-dpdf','-painters',[fig_path expt_ee '_mc_svd_core_nums_' ...
%     ge_type '.pdf'])

%% plot prediction stats
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
% ep
mcs = cell2mat(cos_sim_avg.ep');
scatter((3*stepsz-binsz)*ones(size(mcs(:,1))),mcs(:,1),30,mycc.red_light,'+','linewidth',linew);
scatter((3*stepsz+binsz)*ones(size(mcs(:,2))),mcs(:,2),30,mycc.red,'+','linewidth',linew);
plot([(3*stepsz-binsz)*ones(size(mcs(:,1))),(3*stepsz+binsz)*ones(size(mcs(:,1)))]',...
    mcs','color',mycc.gray);
plot([3*stepsz-binsz*1.5,3*stepsz-binsz*0.5],mean(mcs(:,1))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
plot([3*stepsz+binsz*0.5,3*stepsz+binsz*1.5],mean(mcs(:,2))*ones(2,1),'color',...
    mycc.black,'linewidth',linew);
xlim([0.2 4*stepsz-0.2])
% ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'SVD','CRF','EP'})
ylabel('Similarity')

% accuracy
subplot(1,4,2); hold on
h = boxplot(pred_stats.svd(:,1),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,1),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.ep(:,1),'positions',3*stepsz,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'SVD','CRF','EP'})
ylabel('Accuracy')
box off

% precision
subplot(1,4,3); hold on
h = boxplot(pred_stats.svd(:,2),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,2),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.ep(:,2),'positions',3*stepsz,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'SVD','CRF'})
ylabel('Precision')
box off

% recall
subplot(1,4,4); hold on
h = boxplot(pred_stats.svd(:,3),'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew)
h = boxplot(pred_stats.crf(:,3),'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew)
h = boxplot(pred_stats.ep(:,3),'positions',3*stepsz,'width',ww,'colors',mycc.red);
setBoxStyle(h,linew)
xlim([0 4*stepsz]); ylim([0 1])
set(gca,'xtick',[1,2,3]*stepsz);
set(gca,'xticklabel',{'SVD','CRF','EP'})
ylabel('Recall')
box off

end
