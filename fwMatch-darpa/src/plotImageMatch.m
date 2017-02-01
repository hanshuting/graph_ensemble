function h = plotImageMatch(perm, cornerFileA, imageFileA, cornerFileB, imageFileB)    
    [h, cornerA, cornerB] = plotImageHelper(cornerFileA, imageFileA, cornerFileB, imageFileB);
    
    lineXs = [cornerA(:,1) cornerB(perm,1)];
    lineYs = [cornerA(:,2) cornerB(perm,2)];
           
    % TODO: Map match, proba match.
    
    nSegments = size(lineXs, 1);
    for n = 1:nSegments
        if n == perm(n)
            line(lineXs(n,:), lineYs(n,:), 'Color', 'green');
        else
            line(lineXs(n,:), lineYs(n,:), 'Color', 'red');            
        end
    end

    hold off;        
    

end

