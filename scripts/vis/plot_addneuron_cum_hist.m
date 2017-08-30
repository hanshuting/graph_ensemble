function [] = plot_addneuron_cum_hist(data,mycc,xvec,linew)
% data is a M-by-N cell array, rows contain repeats, columns contain
% different conditions

num_bins = size(data{1,1},1);

if nargin<3
    xvec = 1:num_bins;
end

hold on

plot(xvec,cell2mat(data(:,1)'),'color',mycc.gray,'linewidth',linew);
plot(xvec,cell2mat(data(:,2)'),'color',mycc.orange_light,'linewidth',linew);

h1 = plot(xvec,nanmean(cell2mat(data(:,1)'),2),'color',mycc.black,'linewidth',2*linew);
h2 = plot(xvec,nanmean(cell2mat(data(:,2)'),2),'color',mycc.orange,'linewidth',2*linew);

xlim([xvec(1) xvec(end)]);ylim([-0.01 1.01])
set(gca,'ytick',[0 0.5 1],'linewidth',linew,'fontsize',10)
legend([h1 h2],'original','add neuron'); legend boxoff

end