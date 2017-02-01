% some old code for backup

%% loopy ep subgraph prediction
trilgraph = tril(graph);
num_edge = sum(sum(logical(trilgraph)));
edge_list = zeros(num_edge,2);
[edge_list(:,2),edge_list(:,1)] = find(trilgraph);
G_on = zeros(num_node,num_node);
for i = 1:size(G,2)
    node_1 = edge_list(i,1);
    node_2 = edge_list(i,2);
    G_on(node_1,node_2) = G(4,i);
end

pred_ep_core = zeros(size(data_high,1),1);
for i = 1:length(pred_ep_core)
    fdata = logical(data_high(i,:));
    pred_ep_core(i) = sum(F(2,fdata))+sum(sum(G_on(fdata,fdata)));
end
pred_ep_core_exp = exp(pred_ep_core);
pred_ep_core_nor{e} = pred_ep_core_exp/max(pred_ep_core_exp);

ts = pred_ep_core_nor{e};
ep_pot_pred_thresh = 3*std(ts(ts<quantile(ts,qnoise)));
ep_pot_pred{e} = pred_ep_core_nor{e}>ep_pot_pred_thresh;

% prediction statistics
TP = sum(ep_pot_pred{e}==1&vis_stim_high==ee_stim(e));
TN = sum(ep_pot_pred{e}==0&vis_stim_high~=ee_stim(e));
FP = sum(ep_pot_pred{e}==1&vis_stim_high~=ee_stim(e));
FN = sum(ep_pot_pred{e}==0&vis_stim_high==ee_stim(e));
ee_acc = (TP+TN)/(TP+TN+FN+FP);
ee_prc = TP/(TP+FP);
ee_rec = TP/(TP+FN);
crf_coact_pot_stats{n}(e,:) = [ee_acc,ee_prc,ee_rec];

%% --------- plot all prediction ----------%
    num_frame = length(vis_stim_high);
    num_model = 6;
    vis_stim_mat = repmat(vis_stim_high,num_expt,1));
%     for i = 1:num_expt
%         vis_stim_mat(:,(i-1)*num_frame+1:i*num_frame) = repmat(vis_stim_high'.*...
%             double(vis_stim_high'==i),num_model,1);
%     end
    
    figure;rr = 0.3;
    set(gcf,'color','w','position',[2000,500,800,270],'PaperPositionMode','auto');
%     set(gcf,'PaperSize',[7 8],'PaperPosition',[0 0 7 8]);
    imagesc(vis_stim_mat);
    colormap(cmap);hold on;
    for i = 1:num_expt
        
        % plot raster
        lh = plot((i-1)*num_frame+1+[find(svd_pred{i}),find(svd_pred{i})]',...
            [0.8*ones(sum(svd_pred{i}),1),1*ones(sum(svd_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot((i-1)*num_frame+1+[find(ep_pred{i}),find(ep_pred{i})]',...
            [1.8*ones(sum(ep_pred{i}),1),2*ones(sum(ep_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot((i-1)*num_frame+1+[find(mem_cc_pred{i}),find(mem_cc_pred{i})]',...
            [2.8*ones(sum(mem_cc_pred{i}),1),3*ones(sum(mem_cc_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot((i-1)*num_frame+1+[find(mem_loopy_pred{i}),find(mem_loopy_pred{i})]',...
            [3.8*ones(sum(mem_loopy_pred{i}),1),4*ones(sum(mem_loopy_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot((i-1)*num_frame+1+[find(cent_cc_pred{i}),find(cent_cc_pred{i})]',...
            [4.8*ones(sum(cent_cc_pred{i}),1),5*ones(sum(cent_cc_pred{i}),1)]',...
            'color',[0,0,0,rr]);
        plot((i-1)*num_frame+1+[find(cent_loopy_pred{i}),find(cent_loopy_pred{i})]',...
            [5.8*ones(sum(cent_loopy_pred{i}),1),6*ones(sum(cent_loopy_pred{i}),1)]',...
            'color',[0,0,0,rr]);
%         plot((i-1)*num_frame+1+[find(ep_pot_pred{i}),find(ep_pot_pred{i})]',...
%             [6.5*ones(sum(ep_pot_pred{i}),1),7.5*ones(sum(ep_pot_pred{i}),1)]',...
%             'color',0.7*[1,1,1]);
        
        % plot sim
        plot((i-1)*num_frame+1:i*num_frame,1.4-sim_svd{i}*0.4,'color',0*[1,1,1]);
        plot((i-1)*num_frame+1:i*num_frame,2.4-sim_edge_loopy{i}*0.4,'color',0*[1,1,1])
        plot((i-1)*num_frame+1:i*num_frame,3.4-sim_mem_cc{i}*0.4,'color',0*[1,1,1])
        plot((i-1)*num_frame+1:i*num_frame,4.4-sim_mem_loopy{i}*0.4,'color',0*[1,1,1])
        plot((i-1)*num_frame+1:i*num_frame,5.4-sim_cent_cc{i}*0.4,'color',0*[1,1,1])
        plot((i-1)*num_frame+1:i*num_frame,6.4-sim_cent_loopy{i}*0.4,'color',0*[1,1,1])
%         plot((i-1)*num_frame+1:i*num_frame,7.5-pred_ep_core_nor{i},'color',0*[1,1,1])
        
    end
    ylim([0.8 6.5])
    set(gca,'xtick',[],'ytick',[])
    
    saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.fig']);
    saveas(gcf,[fig_path expt_name{n} '_k_' num2str(k) '_all_prediction.pdf']);

    
         %% --------- plot prediction ----------%
        figure;set(gcf,'color','w','position',[2000,50,930,1050]);
        set(gcf,'PaperSize',[7 8],'PaperPosition',[0 0 7 8]);
        % svd
        subplot(7,1,1);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_svd+0.5,'k','linewidth',1)
        plot([0,length(sim_svd)],1.5-[sim_svd_thresh,sim_svd_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('svd')
        % cc mem
        subplot(7,1,2);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_mem_cc+0.5,'k','linewidth',1)
        plot([0,length(sim_mem_cc)],1.5-[sim_mem_cc_thresh,...
            sim_mem_cc_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('cc mem')
        % loopy mem
        subplot(7,1,3);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_mem_loopy+0.5,'k','linewidth',1)
        plot([0,length(sim_mem_loopy)],1.5-[sim_mem_loopy_thresh,...
            sim_mem_loopy_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('loopy mem')
        % cc cent
        subplot(7,1,4);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_cent_cc+0.5,'k','linewidth',1)
        plot([0,length(sim_cent_cc)],1.5-[sim_cent_cc_thresh,...
            sim_cent_cc_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('cc cent')
        % loopy cent
        subplot(7,1,5);
        imagesc(vis_stim_high')
        hold on;plot(sim_cent_loopy+0.5,'k','linewidth',1)
        plot([0,length(sim_cent_loopy)],1.5-[sim_cent_loopy_thresh,...
            sim_cent_loopy_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('loopy cent')
        % loopy ep cosine
        subplot(7,1,6);
        hold off;imagesc(vis_stim_high')
        hold on;plot(sim_edge_loopy+0.5,'k','linewidth',1)
        plot([0,length(sim_edge_loopy)],1.5-[sim_edge_loopy_thresh,...
            sim_edge_loopy_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('crf ep')
        % loopy ep subgraph
        subplot(7,1,7)
        hold off;imagesc(vis_stim_high')
        hold on;plot(1-pred_ep_core_nor+0.5,'k','linewidth',1)
%         plot([0,length(sim_svd)],1.5-[sim_svd_thresh,sim_svd_thresh],'b--');
        colormap(cmap);set(gca,'yticklabel',{'1' '0.5' '0'});
        title('crf ep potential')
        
        saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core_prediction.fig']);
        saveas(gcf,[fig_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k)...
            '_core_prediction.pdf']);
        
        