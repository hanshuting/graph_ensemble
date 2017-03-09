function [auc,xx_out,yy_out] = plotROC(labels,scores,numClass,cc,linew)
% plot ROC class for multiclass classifiers on current axis handle; return
% AUC vector

if nargin < 4
    % cc = jet(numClass+1);
    cc = jet(numClass);
    cc = max(cc-0.3,0);
end
if nargin < 5
    linew = 1;
end

auc = zeros(numClass,1);
xx_out = [];
yy_out = [];

h = cell(numClass,1);
hstr = '';
lstr = '';
hold on
% patch([0 1 1 0 0],[0 0 1 1 0],0.9*[1 1 1],'edgecolor','none')
plot([0 1],[0 1],'--','color','k','linewidth',1);
for n = 1:numClass
    if sum(labels==n)~=0
        [xx,yy,~,auc(n)] = perfcurve(labels,scores(:,n),n);
        h{n} = plot(xx,yy,'color',cc(n,:),'linewidth',linew);
        hstr = [hstr sprintf('h{%u}',n) ' '];
        lstr = [lstr '''' num2str(n) ''','];
%     text(10,10,sprintf('AUC=%1.2f',auc(n)));
        xx_out(end+1,:) = reshape(xx,1,[]);
        yy_out(end+1,:) = reshape(yy,1,[]);
    end
end
set(gca,'linewidth',1)
xlabel('FPR');ylabel('TPR')

eval(sprintf('legend([%s],%s)\n',hstr(1:end-1),lstr(1:end-1)));

end