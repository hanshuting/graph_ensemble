function cost = matchingCost(C, perm)
% matchingCost
    
    cost = 0;
    D = size(C, 1);
    for i = 1:D
        cost = cost + C(i, perm(i));
    end

end

