function x = solveMAPMosek(gradNodes, gradEdges, Aeq, beq)

%     opts = mskoptimset('Diagnostics','on');

    prob.c = vertcat(gradNodes(:), gradEdges(:));
    prob.a = Aeq;
    
    WIGGLE = 0;
    
    prob.buc = beq + WIGGLE;
    prob.blc = beq - WIGGLE;
    
    prob.blx = [];
    prob.bux = [];

    [r, res] = mosekopt('minimize', prob);
    x = res.sol.bas.xx;

end

