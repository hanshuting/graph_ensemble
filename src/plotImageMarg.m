function h = plotImageMarg(leftNode, rightBeliefs, cornerFileA, imageFileA, cornerFileB, imageFileB)    
% h = plotImageMarg
%
%   Since our Ys are column-permutation matrices, the kth column of beliefs
%   is a distribution for the right-hand-side nodes that the left-hand node
%   k can map to.

    [h, cornerA, cornerB] = plotImageHelper(cornerFileA, imageFileA, cornerFileB, imageFileB);
    
    nNodes = size(cornerA, 1);
    leftXs = repmat(cornerA(leftNode,1), nNodes, 1);
    leftYs = repmat(cornerA(leftNode,2), nNodes, 1);
    
    lineXs = [leftXs cornerB(:,1)];
    lineYs = [leftYs cornerB(:,2)];

    hc = colorbar;
    set(hc, 'YTick', [0 1]);
    set(hc, 'YTickLabel', {''});
    colors = cmapping(rightBeliefs, 'jet');
        
    for n = 1:nNodes
        line(lineXs(n,:), lineYs(n,:), 'Color', colors(n,:));
    end
%     axes('ColorOrder', colors);

    % TODO: COLORBAR SHOULD BE BETWEEN 0 AND 1... FAIL.

    hold off;    

end

