function [] = plot_core_ens(results,coords)
% plot_core_ens(results)
% plot detected ensembles and core neurons

%% parse inputs
auc = results.auc;
auc_ctrl = results.auc_ctrl;
auc_thr = results.auc_thr;
best_model = results.best_model;
core_crf = results.core_crf;
ens_crf = results.ens_crf;
ens_crf_sig = results.ens_crf_sig;
ns_ens = results.ns_ens;
ns_ens_thr = results.ns_ens_thr;
ns_stim = results.ns_stim;
ns_stim_thr = results.ns_stim_thr;
sig_indx = results.sig_indx;

% get parameters
[num_node,num_stim] = size(auc);

%% make coordinates for stimulus nodes
if num_stim > 8
    error('this plotting function cannot make stimulus node coordinates for more than 8 stimli')
else
    mgsz = 10;
    xma = max(coords(:,1))+mgsz;
    yma = max(coords(:,2))+mgsz;
    coords_ext = coords;
    coords_ext(num_node-num_stim+1,:) = [-mgsz -mgsz];
    coords_ext(num_node-num_stim+2,:) = [-mgsz yma];
    coords_ext(num_node-num_stim+3,:) = [xma yma];
    coords_ext(num_node-num_stim+4,:) = [xma -mgsz];
    coords_ext(num_node-num_stim+5,:) = [-mgsz yma/2];
    coords_ext(num_node-num_stim+6,:) = [xma/2 yma];
    coords_ext(num_node-num_stim+7,:) = [xma yma/2];
    coords_ext(num_node-num_stim+8,:) = [xma/2 -mgsz];
    coords_ext = coords_ext(1:num_node,:);
end

% flip coordinates
coords_ext = [coords_ext(:,1),-coords_ext(:,2)];

%% make plots for significant ensemble
mksz = 30;
figure; set(gcf,'color','w')
for n = 1:num_stim

    % 1. plot ensemble with graph structure
    subplot(2,num_stim,n); hold on

    % plot edges
    for ii = 1:num_node
        % set connecting edge to stimulus node red
        if ii == num_node-num_stim+n
            E = find(best_model.graph(ii,:));
            if ~isempty(E)
                linkNode = coords_ext(E,:);
                crNode = repmat(coords_ext(ii,:),length(E),1);
                for jj = 1:length(E)
                    if any(ens_crf_sig{n}==E(jj))
                        plot([crNode(jj,1),linkNode(jj,1)]',[crNode(jj,2),linkNode(jj,2)]',...
                            'color',[1 0 0]);
                    else
                        plot([crNode(jj,1),linkNode(jj,1)]',[crNode(jj,2),linkNode(jj,2)]',...
                            'color',[1 0.7 0.7]);
                    end
                end
            end
        else
            cc = 0.7*[1 1 1];
            E = find(best_model.graph(ii,:));
            if ~isempty(E)
                linkNode = coords_ext(best_model.graph(ii,:)~=0,:);
                crNode = repmat(coords_ext(ii,:),length(E),1);
                for jj = 1:length(E)
                    plot([crNode(jj,1),linkNode(jj,1)]',[crNode(jj,2),linkNode(jj,2)]',...
                        'color',cc);
                end
            end
        end
    end

    % plot nodes
    for ii = 1:num_node
        if ii == num_node-num_stim+n % plot added nodes with square
            cc = [1 0 0];
            mk = 's';
        elseif ii > num_node-num_stim
            cc = 0.7*[1 1 1];
            mk = 's';
        elseif any(ens_crf_sig{n}==ii) % color connected nodes
            cc = [1 0 0];
            mk = 'o';
        elseif any(ens_crf{n}==ii)
            cc = [1 0.7 0.7];
            mk = 'o';
        else
            cc = 0.7*[1 1 1];
            mk = 'o';
        end
        scatter(coords_ext(ii,1),coords_ext(ii,2),mksz,mk,'markerfacecolor',cc,...
            'markeredgecolor','k','linewidth',1);
    end
    axis off equal tight
    
    title(['ensemble #' num2str(n)])
    
    
    % 2. plot scatter
    subplot(2,num_stim,n+num_stim); hold on
    plotGraphHighlight(coords,results.ens_crf{n},[1 0.7 0.7]);
    plotGraphHighlight(coords,results.ens_crf_sig{n},[1 0 0]);
%     title(['significant ensemble #' num2str(n)])
    
end

%% make plots for core ensemble
mksz = 30;
figure; set(gcf,'color','w')
for n = 1:num_stim

    % 1. scatter plots
    subplot(2,num_stim,n); hold on
    
    ns_ens_nor = ns_ens{n}(sig_indx{n})./ns_ens_thr{n}(sig_indx{n});
    auc_nor = auc(ens_crf{n}(sig_indx{n}),n)./auc_thr{n}(ens_crf{n}(sig_indx{n}));
    indx = find((ns_ens_nor>1) & (auc_nor>1));
    nsmi = min(ns_ens_nor);
    nsma = max(ns_ens_nor);
    aucmi = min(auc_nor);
    aucma = max(auc_nor);
    
    % plot threshold
    if ~isempty(aucmi)
        plot([1 1],[aucmi aucma],'k--');
        plot([nsmi nsma],[1 1],'k--');
    end
    
    % plot nodes
    scatter(ns_ens_nor,auc_nor,mksz,[1 0.7 0.7],'o','filled')
    scatter(ns_ens_nor(indx),auc_nor(indx),mksz,[1 0 0],'o','filled')
    
    if ~isempty(aucmi)
        xlim([nsmi nsma]); ylim([aucmi aucma])
    end
    xlabel('normalized node strength')
    ylabel('normalized AUC')
    title(['core #' num2str(n)])
    
    
    % 2. plot cores
    subplot(2,num_stim,n+num_stim); hold on
    plotGraphHighlight(coords,results.ens_crf_sig{n},[1 0.7 0.7]);
    plotGraphHighlight(coords,results.core_crf{n},[1 0 0]);
    
end

end