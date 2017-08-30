function [] = plot_box_seq_single(data,xseq,wr,step,linew,cc)
% for fig6 - plot ensemble reduction add neuron performance

hold on;
for ii = 1:length(xseq)
    h = boxplot(mean(squeeze(data(:,ii,:)),2),'positions',...
        xseq(ii),'width',step*wr,'colors',cc);
    setBoxStyle(h,linew);
end
xlim([xseq(1)-step,xseq(end)+step])
ylim([0 max(data(:))])
set(gca,'xtick',xseq,'xticklabel',xseq,'XTickLabelRotation',45)


end