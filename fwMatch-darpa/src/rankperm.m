function r = rankperm(p)
% rankperm  Return lexicographic rank of permutation

    D = length(p);
    factoraid = zeros(1, D);
    id = 1:D;
    
    % Begin with identity permutation. Iteratively find the indices of p(i)
    % in the identity permutation and remove them.
    
    for i = 1:D
        ip = find(p(i) == id);
        factoraid(i) = ip;
        id(ip) = [];                
    end
    
    factoraid = factoraid - 1; % convert to 0 index
    
    r = 1; % 1 index
    for i = 1:D
        % Position i in the array is radix place value (D - i)!
        place = factorial(D - i);
        r = r + factoraid(i) * place;
    end


end

