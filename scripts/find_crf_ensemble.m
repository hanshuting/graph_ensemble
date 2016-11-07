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
    svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
    num_stim = length(unique(vis_stim))-1;
    num_node = size(best_model.graph,1);
    num_frame = length(Pks_Frame);
    vis_stim_high = vis_stim(Pks_Frame);
    data_high = Spikes(:,Pks_Frame);
    
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
    perc_vec = 0:0.05:0.5;
    acc_vec = zeros(length(perc_vec),num_stim);
    for ii = 1:num_stim
        true_label = double(vis_stim_high==ii)';
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
    
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat'],...
        'core_crf','core_svd');
    
end

end