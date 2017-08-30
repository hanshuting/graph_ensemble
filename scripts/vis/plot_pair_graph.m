function [pval] = plot_pair_graph(vec1,vec2,cc1,cc2,p)
% plot paired data, with linked line and mean per group
% use Wilcoxon signed-rank test for significance test

stepsz = 0.5;
binsz = 0.1;
markersz = 20;
linew = 0.5;
wdr = 0.3;

vec1 = reshape(vec1,[],1);
vec2 = reshape(vec2,[],1);

% [~,pval] = ttest(vec1,vec2);
pval = signrank(vec1,vec2);

hold on
plot([(stepsz-binsz)*ones(size(vec1)),(stepsz+binsz)*...
    ones(size(vec2))]',[vec1,vec2]','color',0.7*[1 1 1]);

scatter((stepsz-binsz)*ones(size(vec1)),vec1,markersz,cc1,'o','filled','linewidth',linew);
scatter((stepsz+binsz)*ones(size(vec2)),vec2,markersz,cc2,'o','filled','linewidth',linew);

plot(stepsz-binsz*[1+wdr 1-wdr],nanmean(vec1)*ones(2,1),'color',cc1,'linewidth',3*linew);
plot(stepsz+binsz*[1+wdr 1-wdr],nanmean(vec2)*ones(2,1),'color',cc2,'linewidth',3*linew);

if pval<=p
    scatter(stepsz,max([vec1;vec2]),'k*');
end

xlim([stepsz-binsz*2 stepsz+binsz*2])

end