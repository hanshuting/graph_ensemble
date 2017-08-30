function [] = crf_find_ensembles(param)
% identify ensembles by switching on and off each neuron
% This code works for multiple stimuli

% parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.ens;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
rwbmap = param.rwbmap;
num_expt = length(expt_name);

load(ccode_path);
load(rwbmap);

%%
for n = 1:num_expt
    
    expt_ee = ee{n}{1};
    
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    
    num_stim = length(unique(setdiff(vis_stim,0)));
    num_node = size(best_model.graph,1);
    num_frame = length(Pks_Frame);
    vis_stim_high = vis_stim(Pks_Frame);
    load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee '.mat']);
    data_high = data';
    
    best_model.ep_on = getOnEdgePot(best_model.graph,best_model.G);
    best_model.ep_on = best_model.ep_on - tril(best_model.ep_on);
    epsum = sum(best_model.ep_on,2);
    epsum(sum(best_model.graph,2)==0) = NaN;
    
    % shuffled models
    shuffle_model = load([model_path 'shuffled_' expt_name{n} '_' ...
        expt_ee '_loopy_fulldata.mat']);
    for ii = 1:length(shuffle_model.graphs)
        shuffle_model.ep_on{ii} = getOnEdgePot(shuffle_model.graphs{ii},...
            shuffle_model.G{ii})';
        shuffle_model.epsum{ii} = sum(shuffle_model.ep_on{ii},2);
        shuffle_model.epsum{ii}(sum(shuffle_model.graphs{ii},2)==0) = NaN;
    end
    shuffle_model.mepsum = nanmean(cellfun(@(x) nanmean(x),shuffle_model.epsum));
    shuffle_model.sdepsum = nanstd(cellfun(@(x) nanmean(x),shuffle_model.epsum));
    
    %% find ensemble with CRF
    % predict each neuron in turn
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
    
    % ------------ AUC - should work for multiple stimuli ------------ %
    auc = zeros(num_node,num_stim);
    for ii = 1:num_stim
        true_label = double(vis_stim_high==ii)';
        for jj = 1:num_node
            [~,~,~,auc(jj,ii)] = perfcurve(true_label,LL_on(jj,:),1);
        end
    end

    auc_ens = cell(num_stim,1);
    core_crf = cell(num_stim,1);
    for ii = 1:num_stim
        num_ens = sum(best_model.graph(num_node-num_stim+ii,:));
        for jj = 1:100
            rd_ens = randperm(num_node,num_ens);
            [~,sim_core] = core_cos_sim(rd_ens,data_high',...
                true_label);
            [~,~,~,auc_ens{ii}(jj)] = perfcurve(true_label,sim_core,1);
        end
        core_crf{ii} = find(auc(:,ii)>(mean(auc_ens{ii})+std(auc_ens{ii}))&...
            (epsum>(shuffle_model.mepsum+shuffle_model.sdepsum)));
        core_crf{ii} = setdiff(core_crf{ii},num_node-num_stim+ii);
    end
    
    
    % --------------------- node strength + AUC --------------------- %
    nodesz = 30;
    nsmi = min(epsum);
    nsma = max(epsum);
    aucmi = 0;
    aucma = 1;
    figure; set(gcf,'color','w','position',[1967 615 555 253])
    for ii = 1:num_stim
        subplot(1,num_stim,ii); hold on
        scatter(epsum,auc(:,ii),nodesz,mycc.gray,'filled')
        scatter(epsum(core_crf{ii}),auc(core_crf{ii},ii),nodesz,mycc.red,'filled')
%         plot([nsmi nsma],th(ii)*[1 1],'k--');
        plot([nsmi nsma],mean(auc_ens{ii})*[1 1],'k--');
        plot([nsmi nsma],(mean(auc_ens{ii})+std(auc_ens{ii}))*[1 1],'--',...
            'color',mycc.gray_light);
        plot([nsmi nsma],(mean(auc_ens{ii})-std(auc_ens{ii}))*[1 1],'--',...
            'color',mycc.gray_light);
        plot(shuffle_model.mepsum*[1 1],[aucmi aucma],'k--');
        plot((shuffle_model.mepsum+shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
            'color',mycc.gray_light);
        plot((shuffle_model.mepsum-shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
            'color',mycc.gray_light);
        xlim([nsmi nsma]); ylim([aucmi aucma])
        xlabel('node strength'); ylabel(['AUC' num2str(ii)]);
    end
    print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_core_NS_AUC.pdf']);
    
    %% plot ensemble highlight
    figure; set(gcf,'color','w','position',[1964 220 864 326])
    for ii = 1:length(core_crf)
        subplot(1,length(core_crf),ii)
        plotGraphHighlight(Coord_active,core_crf{ii},mycc.red)
    end
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
        expt_ee '_vis_core.pdf'])
    
    %% save result
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_core.mat'],'core_crf');
    
end


end

