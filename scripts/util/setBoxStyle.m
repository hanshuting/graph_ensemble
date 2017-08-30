function [] = setBoxStyle(h,linew)

set(h,'linewidth',linew);
set(h(7,:),'visible','off');
set(h(6,:),'linewidth',2*linew);
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')

end