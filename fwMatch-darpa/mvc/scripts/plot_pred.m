function [] = plot_pred(pred,sim,sim_thresh,stim_vec,cmap)

num_frame = length(stim_vec);
num_model = size(pred,1);
vis_stim_mat = repmat(stim_vec,1,num_model)';

rr = 0.5;
ss = 0.1;

figure;
set(gcf,'color','w')
set(gcf,'position',[2013 80 1171 671]);
set(gcf,'PaperPositionMode','auto');

imagesc(vis_stim_mat);
colormap(cmap);
hold on;

for n = 1:num_model

    % plot raster
    plot([find(pred(n,:));find(pred(n,:))],[(0.5*(n+1)+ss)*ones(sum(pred(n,:)),1),...
        (0.5*(n+1)-ss)*ones(sum(pred(n,:)),1)]','color',[0,0,0,rr]);
    
    % plot sim
    plot(1:num_frame,0.5*(n+1)+num_model*0.5-0.5*sim(n,:),'color',0*[1,1,1]);
    plot([1,num_frame],0.5*(n+1)+num_model*0.5-0.5*sim_thresh(n)*[1,1],'-','color',0.8*[1,1,1]);
    
end
ylim([1-ss num_model+0.5])
set(gca,'tickdir','out','ytick',[])


end