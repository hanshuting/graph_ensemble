function [] = plot_opto_spont_cum_hist(data,mycc,xvec,linew)
% data is a M-by-N cell array, rows contain repeats, columns contain
% different conditions

num_bins = size(data{1,1},1);

if nargin<3
    xvec = 1:num_bins;
end

hold on

plot(xvec,log(cell2mat(data(:,1)')),'color',mycc.gray,'linewidth',linew);
plot(xvec,log(cell2mat(data(:,2)')),'color',mycc.blue_light,'linewidth',linew);

h1 = plot(xvec,log(nanmean(cell2mat(data(:,1)'),2)),'color',mycc.black,'linewidth',2*linew);
h2 = plot(xvec,log(nanmean(cell2mat(data(:,2)'),2)),'color',mycc.blue,'linewidth',2*linew);

ymin = log(min(reshape(cell2mat(data),[],1)));
xlim([xvec(1) xvec(end)]);ylim([ymin 0])
set(gca,'linewidth',linew,'fontsize',10)
legend([h1 h2],'pre','post'); legend boxoff

end