function diffs = featureTest(M, survey)

    E = size(M, 1);
    diffs = zeros(E, 1);
    for e = 1:E
        xi = survey(M(e,1),:);
        xj = survey(M(e,2),:);
        diffs(e) = norm(xi - xj, 1);
    end
    
end

