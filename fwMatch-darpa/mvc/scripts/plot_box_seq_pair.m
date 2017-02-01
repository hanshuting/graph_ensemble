function [] = plot_box_seq_pair(data,xseq,binsz,wr,step,linew,mycc)
% for fig6 - plot ensemble reduction add neuron performance

hold on;
for ii = 1:length(xseq)
    h = boxplot(mean(squeeze(data(:,ii,:,1)),2),'positions',...
        xseq(ii)-binsz,'width',step*wr,'colors',mycc.gray);
    setBoxStyle(h,linew);
    h = boxplot(mean(squeeze(data(:,ii,:,2)),2),'positions',...
        xseq(ii)+binsz,'width',step*wr,'colors',mycc.orange);
    setBoxStyle(h,linew);
end
xlim([xseq(1)-step,xseq(end)+step])
ylim([0 max(data(:))])
set(gca,'xtick',xseq,'xticklabel',xseq,'XTickLabelRotation',45)

end