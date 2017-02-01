function [ h ] = paperFig(wf, hf)
% h = paperFig(wf, hf) create publication figure of width wf and hf
% relative to _text width_ with a 1 in magin on 8.5 x 11in page.
    margin = 1;
    pageWidth = 8.5;
    textWidth = pageWidth - 2*margin;
    
    w = wf * textWidth;
    h = hf * textWidth;
    
    h = fig('units','inches','width',w,'height',h,'font','Helvetica','fontsize',10);
    set(h,'PaperPositionMode','auto')
end

