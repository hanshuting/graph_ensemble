function [logZ, dLogZ] = exactLogZNoFeat(A)

    D = size(A, 1);
    Z = ryser_small(exp(A));
    logZ = log(Z);
    
    % d \log Z / d A_{ij} = E[\exp(A_{ij})] (log-linear model)    
    dLogZ = factorial(D - 1) / Z * exp(A);
    
    Ys = perms(1:D);    
    M = size(Ys, 1);
    
    % Sum over each permutation.    
    dLogZ = zeros(D, D);
    Z = 0;
    for m = 1:M
        Ym = expandPerm(Ys(m,:));                
        expTrAY = exp(sum(vec(A(Ym))));
        Z = Z + expTrAY;
        for j = 1:D
            i = Ys(m,j);
            dLogZ(i,j) = dLogZ(i,j) + expTrAY;
        end        
    end
        
    logZ = log(Z);
    dLogZ = dLogZ / Z;
end

