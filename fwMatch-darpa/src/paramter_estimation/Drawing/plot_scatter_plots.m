function plot_scatter_plots(thetaN1,thetaN2,thetaE1,thetaE2,fid,ftitle,edges)

    a  =40; %radius of the points
    clf;
    figure(fid),
    title(ftitle),
    subplot(1,2,1),
    hold on,
    xlim([min(thetaN2)-0.00001 max(thetaN2)+0.00001]),
    ylim([min(thetaN1) max(thetaN1)]),
    ptN1 = linspace(min(thetaN1),max(thetaN1),50);
    ptN2 = linspace(min(thetaN2),max(thetaN2),50);
    % Y = x
    plot(ptN1,ptN1,'Color','r','LineWidth',2);
    % x = 0
%     plot(ptN2,zeros(1,length(ptN2)),'Color','b','LineWidth',2);
%     % y = 0
%     plot(zeros(1,length(ptN1)),ptN1,'Color','g','LineWidth',2);
    scatter(thetaN2,thetaN1,a,'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5),
    for i=1:length(thetaN1)
        text(thetaN2(i),thetaN1(i),['   ' num2str(i)]);
    end
    title('Node Potentials'),
    xlabel('{\theta_N}');
    ylabel('{\theta_N}^*');
    
    subplot(1,2,2),
    hold on,
    xlim([min(thetaE2)-0.00001 max(thetaE2)+0.00001]),
    ylim([min(thetaE1) max(thetaE1)]),
    ptE1 = linspace(min(thetaE1),max(thetaE1),50);
    ptE2 = linspace(min(thetaE2),max(thetaE2),50);
%     plot(ptE1,ptE1,'Color','r','LineWidth',2);
%     plot(ptE2,zeros(1,length(ptE2)),'Color','b','LineWidth',2);
%     plot(zeros(1,length(ptE1)),ptE1,'Color','g','LineWidth',2);
    scatter(thetaE2,thetaE1,a,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
    for i=1:length(thetaE1)
        text(thetaE2(i),thetaE1(i),['   ( ' num2str(edges(i,1)) ', ' num2str(edges(i,2)) ' )']);
    end
    title('Edge Potentials'),
    xlabel('{\theta_E}');
    ylabel('{\theta_E}^*');
    
end

