function [] = plot_pred(pred,sim,sim_thresh,stim_vec,cmap)
% plot_pred(pred,sim,sim_thresh,stim_vec,cmap)
% plot prediction values, threshold, and prediction raster
% INPUT:
%     pred: T-by-num_model binary matrix
%     sim: T-by-num_model prediction values
%     sim_thresh: 1-by-num_model prediction threshold
%     stim_vec: T-by-num_model ground truth labels
%     cmap: color map, used for coloring stimulus stripes in the background
% 

[num_frame, num_model] = size(pred);
stim_mat = repmat(stim_vec,1,num_model)';

rr = 0.5; % transparency
ss = 0.1; % raster height

figure;
set(gcf,'color','w')
% set(gcf,'position',[2013 80 1171 671]);
set(gcf,'PaperPositionMode','auto');

imagesc(stim_mat);
colormap(cmap);
hold on;

for n = 1:num_model

    % plot raster
    pred_time = find(pred(:,n));
    plot([1;1]*pred_time',(0.5*(n+1)+ss*[1;-1])*ones(size(pred_time')),'color',[0,0,0,rr]);
    
    % plot sim
    plot(1:num_frame,0.5*(n+1)+num_model*0.5-0.5*sim(:,n),'color',0*[1,1,1]);
    plot([1,num_frame],0.5*(n+1)+num_model*0.5-0.5*sim_thresh(n)*[1,1],'-','color',0.8*[1,1,1]);
    
end
ylim([1-ss num_model+0.5])
set(gca,'tickdir','out','ytick',[])
xlabel('frame')

end