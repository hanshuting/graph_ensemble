function [] = print_graph_to_pdf(graph, filename)
    % hack to print a biograph to pdf (see http://www.mathworks.com/matlabcentral/newsreader/view_thread/169060) 
    gui = biograph.bggui(graph);
    fig = figure();
    set(fig,'PaperPositionMode','manual');         
    set(fig,'PaperUnits','inches');
    set(fig,'PaperPosition',[0.5,0.5,10,7.5]);
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperType','usletter');
    copyobj(gui.biograph.hgAxes,fig);
    print(fig, '-dpdf', filename);

