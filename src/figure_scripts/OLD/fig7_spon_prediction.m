% script for figure 4, identifying ensembles using different methods

%% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
vis_ee = {{'01_high','02_high'}};
spon_ee = {{'spon_pre','spon_post'}};
num_shuff = 100;
k = 3;
p = [0.05,0.01];
qnoise = 0.5;
save_path = 'C:\Shuting\fwMatch\results\stats\';
all_fig_path = 'C:\Shuting\fwMatch\results\fig\stats\';

%% ensemble predictions
vis_stim_seq = {[1,2]};
num_ee = length(expt_name);

sim_spon_all = cell(num_ee,2,2); % expt-by-vis-by-spon
pred_spon_all = cell(num_ee,2,2);
ov_core_all = cell(num_ee,2,2);
ov_noncore_all = cell(num_ee,2,2);

for n = 1:length(expt_name)
    
    expt_ee = vis_ee{n};
    expt_spon_ee = spon_ee{n};
    num_expt = length(expt_ee);
    load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} '.mat']);
    load(['C:\Shuting\fwMatch\data\ensembles\' expt_name{n} '_core_svd.mat']);
    num_node = size(Spikes,1);
    
    % initialize
    sim_spon = cell(num_expt,2); % store pre in column 1, post in column 2
    pred_spon = cell(num_expt,2);
    sim_spon_thresh = zeros(num_expt,2);
    ov_core = cell(num_expt,2);
    ov_noncore = cell(num_expt,2);
    
    ee_stim = vis_stim_seq{n};
    
    %% ----------------- prediction ---------------- %
    for e = 1:num_expt
        
        % set path
        data_path = ['C:\Shuting\fwMatch\data\' expt_name{n} '\']; 
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\']; 
        cc_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\cc\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        fig_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\fig\'];
        
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
        load([model_path  expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\Pks_Frames.mat']);
        load(['C:\Shuting\fwMatch\data\' expt_name{n} '\' expt_name{n} ...
            '_' expt_ee{e} '.mat']);
        load('C:\Shuting\fwMatch\results\wrbmap.mat');
        
        spon_data = cell(2,1);
        
        for i = 1:2
            
            if ~isempty(expt_spon_ee{i})
                
                load([data_path expt_name{n} '_' expt_spon_ee{i} '.mat']);
                spon_data{i} = data;
                
                % cosine similarity
                core_vec = zeros(num_node,1);
                core_vec(core_edge_loopy) = 1;
                sim_spon{e,i} = 1-pdist2(data,core_vec','cosine');
                ts = sim_spon{e,i};
                sim_spon_thresh(e,i) = 3*std(ts(ts<=quantile(ts,qnoise)))+...
                    mean(ts(ts<=quantile(ts,qnoise)));
                pred_spon{e,i} = ts>sim_spon_thresh(e,i);
                fpred = data(pred_spon{e,i},:);
                ov_core{e,i} = sum(fpred(:,logical(core_vec)),1);
                ov_noncore{e,i} = sum(fpred(:,~logical(core_vec)),1);
                
            else
                
                sim_spon{e,i} = [];
                sim_spon_thresh(e,i) = [];
                pred_spon{e,i} = [];
                ov_core{e,i} = [];
                ov_noncore{e,i} = [];
                
            end
            
        end
        
    end

    % store results
    sim_spon_all(n,:,:) = sim_spon;
    pred_spon_all(n,:,:) = pred_spon;
    ov_core_all(n,:,:) = ov_core;
    ov_noncore_all(n,:,:) = ov_noncore;
    
    %% --------- plot similarity and prediction raster ----------%
    num_frame = size(spon_data{1},1)+size(spon_data{2},1);
    pred_mat = cell2mat(pred_spon')';
    
    rr = 0.3;ss = 0.1;
    figure;
    set(gcf,'color','w','position',[2000,278,800,160],...
        'PaperPositionMode','auto');
    hold on;
    imagesc(1-pred_mat);colormap(gray);
    plot(size(spon_data{1},1)*[1,1],[0.5,2.5],'--','color','r')
    xlim([1 num_frame]);ylim([0.5 2.5]);
    set(gca,'tickdir','out','ytick',[1,2],'yticklabel',{'vis 01','vis 02'})
    box on
    
%     for i = 1:num_expt
%         plot([find(svd_pred{i}),find(svd_pred{i})]',...
%             [(0.5+ss+(i-1)*num_model/2)*ones(sum(svd_pred{i}),1),...
%             (1-ss+(i-1)*num_model/2)*ones(sum(svd_pred{i}),1)]',...
%             'color',[0,0,0,rr]);
%     end
    
    
%     saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.fig']);
%     saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.pdf']);

end

%% plot ensemble overlay in predicted frames
binsz = 0.2;
figure;
set(gcf,'color','w','position',[2000,441,354,260],'PaperPositionMode','auto')
hold on;
h = boxplot(reshape(cell2mat(ov_core_all(:,1,:)),[],1),...
    'positions',1-binsz,'width',binsz*0.5,'colors','r');
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_noncore_all(:,1,:)),[],1),...
    'positions',1+binsz,'width',binsz*0.5,'colors',[0.5 0.5 0.5]);
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_core_all(:,2,:)),[],1),...
    'positions',2-binsz,'width',binsz*0.5,'colors','b');
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_noncore_all(:,2,:)),[],1),...
    'positions',2+binsz,'width',binsz*0.5,'colors',[0.5 0.5 0.5]);
set(h(7,:),'visible','off')
xlim([0.5 2.5]);
ylim([0 max([max(reshape(cell2mat(ov_core_all(:)'),[],1));...
    max(reshape(cell2mat(ov_noncore_all(:)'),[],1))])])
ylabel('# overlay')
set(gca,'xtick',[1-binsz,1+binsz,2-binsz,2+binsz],'xticklabel',...
    {'hor core','hor noncore','ver core','ver noncore'},...
    'XTickLabelRotation',45)
box off

saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_spon_overlay.fig']);
saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_spon_overlay.pdf']);

%% plot ensemble overlay in before and after vis stim
binsz = 0.2;

figure;
set(gcf,'color','w','position',[2000,441,354,260],'PaperPositionMode','auto')
hold on;
h = boxplot(reshape(cell2mat(ov_core_all(:,:,1)),[],1),...
    'positions',1-binsz,'width',binsz*0.5,'colors','k');
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_noncore_all(:,:,1)),[],1),...
    'positions',1+binsz,'width',binsz*0.5,'colors',0.5*[1 1 1]);
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_core_all(:,:,2)),[],1),...
    'positions',2-binsz,'width',binsz*0.5,'colors',[1 0.5 0]);
set(h(7,:),'visible','off')
h = boxplot(reshape(cell2mat(ov_noncore_all(:,:,2)),[],1),...
    'positions',2+binsz,'width',binsz*0.5,'colors',0.5*[1 1 1]);
set(h(7,:),'visible','off')
xlim([0.5 2.5]);
ylim([0 max([max(reshape(cell2mat(ov_core_all(:)'),[],1));...
    max(reshape(cell2mat(ov_noncore_all(:)'),[],1))])])
ylabel('# overlay')
set(gca,'xtick',[1-binsz,1+binsz,2-binsz,2+binsz],'xticklabel',...
    {'pre core','pre noncore','post core','post noncore'},...
    'XTickLabelRotation',45)
box off

saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_spon_overlay.fig']);
saveas(gcf,[all_fig_path 'k_' num2str(k) '_core_spon_overlay.pdf']);

%% save results
save([save_path '_all_pred_stats.mat'],'svd_stats','crf_coact_stats',...
    'cc_comm_stats','crf_comm_stats','cc_cent_stats','crf_cent_stats',...
    '-v7.3');


