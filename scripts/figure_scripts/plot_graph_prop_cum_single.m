function [] = plot_graph_prop_cum_single(data,mycc,xvec,linew)
% for figure3

num_bins = max(cellfun('length',reshape(struct2cell(data),[],1)));
data_cc = cellfun(@(x) padarray(x,num_bins-length(x),NaN,'post'),...
    getNestedField(data,'cc'),'uniformoutput',false);
data_loopy = cellfun(@(x) padarray(x,num_bins-length(x),NaN,'post'),...
    getNestedField(data,'crf'),'uniformoutput',false);
data_cc_rand = cellfun(@(x) padarray(x,num_bins-length(x),NaN,'post'),...
    getNestedField(data,'cc_rand'),'uniformoutput',false);
data_loopy_rand = cellfun(@(x) padarray(x,num_bins-length(x),NaN,'post'),...
    getNestedField(data,'crf_rand'),'uniformoutput',false);

if nargin<3
    xvec = 1:num_bins;
end

hold on

plot(xvec,log(cell2mat(data_cc')),'color',mycc.green_light,'linewidth',linew);
plot(xvec,log(cell2mat(data_loopy')),'color',mycc.orange_light,'linewidth',linew);
plot(xvec,log(cell2mat(data_cc_rand')),'--','color',mycc.gray_light,'linewidth',linew);
plot(xvec,log(cell2mat(data_loopy_rand')),'color',mycc.gray_light,'linewidth',linew);

h1 = plot(xvec,log(nanmean(cell2mat(data_cc'),2)),'color',mycc.green,'linewidth',2*linew);
h2 = plot(xvec,log(nanmean(cell2mat(data_loopy'),2)),'color',mycc.orange,'linewidth',2*linew);
h3 = plot(xvec,log(nanmean(cell2mat(data_cc_rand'),2)),'--','color',mycc.black,'linewidth',2*linew);
h4 = plot(xvec,log(nanmean(cell2mat(data_loopy_rand'),2)),'color',mycc.black,'linewidth',2*linew);

xlim([xvec(1) xvec(end)]);% ylim([-0.01 1.01]
set(gca,'linewidth',linew,'fontsize',10)
legend([h1 h2 h3 h4],'CC','CRF','CCrand','CRFrand'); legend boxoff

end