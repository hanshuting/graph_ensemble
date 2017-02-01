function [h1, h2] = plotVectorSeries(t1, m1, t2, m2)

    % TODO: Generalize to > 2
    
    
    lb = min([m1(:) ; m2(:)]);
    ub = max([m1(:) ; m2(:)]);
    
    h1 = figure; imagesc(m1); caxis([lb, ub]); colorbar; title(t1); xlabel('Iter');
    h2 = figure; imagesc(m2); caxis([lb, ub]); colorbar; title(t2); xlabel('Iter');


end

