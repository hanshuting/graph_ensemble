function [Aeq, beq] = reduceEqualityConstraints(badAeq, badbeq)

    [N, M] = size(badAeq);
    
    badA = [badAeq badbeq];
    [R, jb] = rref(badA);
    
    r = length(jb);
    Aeq = R(1:r,1:M);
    beq = R(1:r,M+1);
    
    assert(all(vec(R(r+1:end,:)) == 0));

end

