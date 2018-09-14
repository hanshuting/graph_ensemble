function [] = fig_find_louvain_ensemble(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
result_path_base = param.result_path_base;
savestr = param.savestr;

%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    
    stim_type = setdiff(unique(vis_stim),0);
    num_stim = length(stim_type);
    
    %% find louvain communities
    gamma = 1:0.05:1.5; num_gamma = length(gamma);
    
    spike_mat = Spikes(:,Pks_Frame);
    dims = size(spike_mat);
    
    % detect ensembles
    M = cell(num_gamma,1); mu = zeros(num_gamma,1);
    CIJ = double(best_model.graph(1:dims(1)-num_stim,1:dims(1)-num_stim));
    for ii = 1:num_gamma
        [M{ii}, mu(ii)] = clusterConsensusLouvainAdjmat(CIJ,gamma(ii));
    end

    [~,best_indx] = min(mu);
    M_best = M{best_indx};

    %% plot
    ulabel = setdiff(unique(M_best),0);
    num_comm = length(ulabel);
    mm = floor(sqrt(num_comm)); nn = ceil(num_comm/mm);
    figure; set(gcf,'color','w')
    for ii = 1:num_comm
        subplot(mm,nn,ii);
        plotGraphHighlight(Coord_active,find(M_best==ulabel(ii)),[1 0 0], 0.5)
        set(gca,'xtick',[],'ytick',[],'ztick',[])
        title(['ensemble #' num2str(ii)]);
    end
    
    %% find visual ensembles
    roc_core_vis = cell(num_comm,1);
    auc_vis = zeros(num_comm,num_stim);
    for ii = 1:num_comm
        core_vec = zeros(1,dims(1));
        core_vec(M_best==ulabel(ii)) = 1;
        for jj = 1:num_stim
            sim_core = 1-pdist2(spike_mat',core_vec,'cosine')';
            [xx,yy,~,auc_vis(ii,jj)] = perfcurve(double(vis_stim(Pks_Frame)==stim_type(jj)),sim_core,1);
            roc_core_vis{ii}{jj} = [xx,yy];
        end
    end

    louvain_core = cell(num_stim,1);
    for ii = 1:num_stim
        [~,indx] = max(auc_vis(:,ii));
        louvain_core{ii} = find(M_best==ulabel(indx));
    end
    
    %% save
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_louvain_core.mat'],'louvain_core');
    
    
end

end