function discrete = rankperms(perms)
% permHist  Plot histogram of permutations

    [N, D] = size(perms);
    discrete = zeros(N, 1);
    
    for n = 1:N
        discrete(n) = rankperm(perms(n,:));
    end


end

