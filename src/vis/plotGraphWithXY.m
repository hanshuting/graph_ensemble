function [h] = plotGraphWithXY(adjmat,coords,weightmat)
% plot adjmat using the given node coordinates, with color specifying the
% weight (red: positive weight; blue: negative weight)

nodesz = 300;
ccRange = max(abs(weightmat(:)))*1.01;
cmap = jet(100);
N = size(adjmat,1);

% extend coordinates
coords(end+1,:) = [0 max(coords(:,2))];
coords(end+1,:) = [0 0];
coords(end+1,:) = [max(coords(:,1)) 0];
coords(end+1,:) = [max(coords(:,1)) max(coords(:,2))];
coords = coords(1:N,:);

% flip coordinates
coords = [coords(:,1),-coords(:,2)];

h = figure;
set(gcf,'color','w');
hold on;
for i = 1:N
    
    % plot edge
    E = find(adjmat(i,:));
    if ~isempty(E)
        linkNode = coords(adjmat(i,:)~=0,:);
        crNode = repmat(coords(i,:),length(E),1);
%         plot([crNode(:,1),linkNode(:,1)]',[crNode(:,2),linkNode(:,2)]','k--');
        for j = 1:length(E)
            cc = cmap(ceil((weightmat(i,E(j))+ccRange)/(2*ccRange)*100),:);
%            cc = cmap(ceil(weightmat(i,E(j))/ccRange*100),:);
            plot([crNode(j,1),linkNode(j,1)]',[crNode(j,2),linkNode(j,2)]',...
                'color',cc);
        end
    end

end

colorbar;colormap(jet);
caxis([-ccRange ccRange]);

% plot node
for i = 1:N
    
    scatter(coords(i,1),coords(i,2),nodesz,'k','linewidth',2);
    
    % show node name
%     text(coords(i,1)-2,coords(i,2),num2str(i),'color','k','fontsize',10,'fontweight','bold');
    
    % or show node degree
%     text(coords(i,1)-2,coords(i,2),num2str(sum(adjmat(i,:))),'color','k','fontsize',10,'fontweight','bold');
    
    % or show node identity
    text(coords(i,1)-2,coords(i,2),num2str(i),'color','k','fontsize',10,'fontweight','bold');
    
end

axis off

end
