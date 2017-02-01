function [] = plot_pred_raster(pred,stim_vec,cmap)

num_frame = length(stim_vec);
num_model = size(pred,1);
vis_stim_mat = repmat(stim_vec,1,num_model)';

rr = 1;
ss = 0.2;

figure;
set(gcf,'color','w')
set(gcf,'position',[2013 421 787 330]);
set(gcf,'PaperPositionMode','auto');

imagesc(vis_stim_mat);
colormap(cmap);
hold on;

for n = 1:num_model

    % plot raster
%     plot([find(pred(n,:));find(pred(n,:))],[(0.5*(n+1)+ss)*ones(sum(pred(n,:)),1),...
%         (0.5*(n+1)-ss)*ones(sum(pred(n,:)),1)]','color',[0,0,0,rr]);
    
    pdt = find(pred(n,:));
    ypos = 0.5*(n+1);
    for ii = 1:length(pdt)
        patch([pdt(ii)-0.5 pdt(ii)+0.5 pdt(ii)+0.5 pdt(ii)-0.5 pdt(ii)-0.5],...
            [ypos-ss ypos-ss ypos+ss ypos+ss ypos-ss],[0 0 0],'edgecolor',...
            'none','facealpha',rr);
    end
    
end
ylim([1-ss 0.5*(n+1)+ss])
set(gca,'tickdir','out','ytick',[])


end