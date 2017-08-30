function [] = find_crf_ensemble(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
num_expt = length(expt_name);

qnoise = 0.7;

load(ccode_path);
load(rwbmap);

%% find ensembles
for n = 1:num_expt
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(setdiff(unique(vis_stim),0));
    num_node = size(best_model.graph,1);
    num_frame = length(Pks_Frame);
    vis_stim_high = vis_stim(Pks_Frame);
    data_high = Spikes(:,Pks_Frame);
    
    % SVD
    core_svd = cell(num_stim,1);
    if exist([data_path 'ensembles\' expt_name{n} '_core_svd.mat'])
        svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
        for ii = 1:num_stim
            for jj = 1:length(svd_data.svd_state)
                if strcmp(num2str(ii),svd_data.svd_state{jj})
                    core_svd{ii} = svd_data.core_svd{jj};
                    break;
                end
            end
        end
    end

    LL_frame = zeros(num_node,num_frame,2);
    for ii = 1:num_node
        for jj = 1:num_frame
            frame_vec = data_high(:,jj)';
            frame_vec(ii) = 0;
            LL_frame(ii,jj,1) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,frame_vec);
            frame_vec(ii) = 1;
            LL_frame(ii,jj,2) = compute_avg_log_likelihood(best_model.node_pot,...
                best_model.edge_pot,best_model.logZ,frame_vec);
        end
    end
    LL_on = squeeze(LL_frame(:,:,2)-LL_frame(:,:,1));
    
    thr = zeros(num_node,1);
    LL_pred = nan(num_node,num_frame);
    for ii = 1:num_node
        [LL_pred(ii,:),thr(ii)] = pred_from_LL(LL_on(ii,:),qnoise);
    end
    LL_pred_count = zeros(num_node,num_stim+1);
    for ii = 1:num_stim
        LL_pred_count(:,ii) = nansum(LL_pred(:,vis_stim_high==ii),2)/sum(vis_stim_high==ii);
    end
    LL_pred_count(:,end) = nansum(LL_pred(:,vis_stim_high==0),2)/sum(vis_stim_high==0);
    [max_perc,node_idt] = max(LL_pred_count,[],2);
    
    % find percentage threshold
    perc_vec = 0.2:0.05:0.4;
    acc_vec = zeros(length(perc_vec),num_stim);
    for ii = 1:num_stim
        true_label = double(vis_stim_high==ii);
        for jj = 1:length(perc_vec)
            core = find(node_idt.*(max_perc>perc_vec(jj))==ii);
            [~,~,~,~,acc] = core_cos_sim(core,data_high',true_label);
            acc_vec(jj,ii) = acc;
        end
    end
    [~,perc_thresh] = max(mean(acc_vec,2));
    perc_thresh = perc_vec(perc_thresh);
    
    core_crf = cell(num_stim+1,1);
    for ii = 1:num_stim+1
        core_crf{ii} = find(node_idt.*(max_perc>perc_thresh)==ii);
    end
    
    %% plot counts
    figure; plot(perc_vec,acc_vec);
    
    dsz = 20;
    stepsz = 0.5;
    cc_dot = {mycc.red,mycc.red_light; mycc.blue,mycc.blue_light;...
        mycc.green,mycc.green_light; mycc.purple,mycc.purple_light};
    figure; set(gcf,'color','w','position',[2271 387 250 277])
    hold on
    noncore = setdiff(1:num_node,cell2mat(core_crf(1:num_stim)))';
    for ii = 1:num_stim-1 % plot non core gray line
        plot(stepsz*[ii;ii+1]*ones(size(noncore))',LL_pred_count(noncore,ii:ii+1)',...
            'color',mycc.gray_light);
    end
    for ii = 1:num_stim % plot core colored ling
        for jj = 1:num_stim-1
            plot(stepsz*[jj;jj+1]*ones(size(core_crf{ii}))',LL_pred_count(core_crf{ii},jj:jj+1)',...
                'color',cc_dot{ii,2});
        end
    end
    for ii = 1:num_stim % plot non core gray dots
        scatter(ii*stepsz*ones(size(noncore)),LL_pred_count(noncore,ii),...
            dsz,mycc.gray,'filled')
    end
    for ii = 1:num_stim % plot core colored dots
        scatter(ii*stepsz*ones(size(core_crf{ii})),LL_pred_count(core_crf{ii},ii),...
            dsz,cc_dot{ii,1},'filled')
    end
    plot([0.3,num_stim*stepsz+0.2],perc_thresh*ones(1,2),'k--');
    xlim([0.3 num_stim*stepsz+0.2])
    set(gca,'xtick',stepsz*(1:num_stim),'xticklabel',1:num_stim)
    ylabel('prediction (%)')
    
    %% save results
    if exist([result_path_base '\' expt_name{n} '\core\'])~=7
        mkdir([result_path_base '\' expt_name{n} '\core\']);
    end
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat'],...
        'core_crf','core_svd');
    
end

end