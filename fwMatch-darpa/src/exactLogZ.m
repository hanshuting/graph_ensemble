function [logZ, dLogZ] = exactLogZ(features, theta)
% exactZ  Compute exact Z, features.
%
%   [Z, dZ] = exactZ(Fkm, theta)
%    - Fm is collection of feature matrices at sample m

    K = length(features);
    D = size(features{1}, 1);
    
    Ys = perms(1:D);    
    M = size(Ys, 1);
    
    sumThetaK = zeros(M, 1);    
    trFY      = zeros(M, K);
    dZ        = zeros(K, 1);
    
    % Sum over each permutation.
    for m = 1:M
        Ym = expandPerm(Ys(m,:));        

        for k = 1:K            
            trFY(m,k) = sum(vec(features{k}(Ym)));
            sumThetaK(m) = sumThetaK(m) + theta(k) * trFY(m,k);                                    
        end                
    end
            
    logZ = logSumExp(sumThetaK);        
    assert(imag(logZ) == 0, 'logZ was imaginary!');
    
    % Compute dZ (see writeup)    
    dZ = sum(bsxfun(@times, exp(sumThetaK), trFY), 1)';
        
    dLogZ = dZ / exp(logZ); % Chain rule
end

