function [objVal, params] = crfObjVal(YN, YE, TN, TE, Ut, Vt, edges, lambda, rho)

    [nNodes, L] = size(TN);

    assert(size(edges, 1) == 2, 'Edges must be 2 x nEdges');    
    nEdges = size(edges, 2);
    nNodes = size(TN, 1);
    edgeRhos = rho * ones(nEdges, 1);
    [RN, RE] = makeRhoMatOvercomplete(edgeRhos, edges, nNodes, L);

    WN = TN - YN;
    WE = TE - YE;
    
    objVal = 1/(2*lambda) * sum(vec((Ut * WN) .^ 2)) + ...
             1/(2*lambda) * sum(vec((Vt * WE) .^ 2)) + ...
             frobProd(1 - RN, TN .* log(TN)) + ...
             frobProd(RE,     TE .* log(TE));
         
    params.F = -(1/lambda) * Ut * WN;
    params.G = -(1/lambda) * Vt * WE;
    

end

