function [] = plotEnsembleIdentification(results)

auc = results.auc;
auc_ctrl = results.auc_ctrl;
ens_crf = results.ens_crf;
core_crf = results.core_crf;
ns = results.ns;
ns_shuff = results.ns_shuff;
nstim = size(auc,2);

%% make figure
% if nargin < 2
    figure; set(gcf,'color','w')
% else
%     figure(h);
% end

%% process results
% normalize auc
auc_nor = auc./mean(auc_ctrl,3);
auc_thr = 1;

% normalize ns
ns_nor = ns./mean(ns_shuff,2);
ns_thr = 1;

%% plot scatter
nodesz = 30;
nsmi = min(ns_nor);
nsma = max(ns_nor);
for ii = 1:nstim
    subplot(1,nstim,ii); hold on
    aucmi = min(auc_nor(:,ii));
    aucma = max(auc_nor(:,ii));
    scatter(ns_nor,auc_nor(:,ii),nodesz,0.7*[1 1 1],'filled')
    scatter(ns_nor(ens_crf{ii}),auc_nor(ens_crf{ii},ii),nodesz,[1 0.8 0.8],'filled')
    scatter(ns_nor(core_crf{ii}),auc_nor(core_crf{ii},ii),nodesz,[1 0.2 0.2],'filled')
    plot([nsmi nsma],auc_thr*[1 1],'k--');
%     plot([nsmi nsma],(mean(auc_ens{ii})+std(auc_ens{ii}))*[1 1],'--',...
%         'color',mycc.gray_light);
%     plot([nsmi nsma],(mean(auc_ens{ii})-std(auc_ens{ii}))*[1 1],'--',...
%         'color',mycc.gray_light);
    plot(ns_thr*[1 1],[aucmi aucma],'k--');
%     plot((shuffle_model.mepsum+shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
%         'color',mycc.gray_light);
%     plot((shuffle_model.mepsum-shuffle_model.sdepsum)*[1 1],[aucmi aucma],'--',...
%         'color',mycc.gray_light);
    xlim([nsmi nsma]); ylim([aucmi aucma])
    xlabel('node strength'); ylabel(['AUC' num2str(ii)]);
end

end