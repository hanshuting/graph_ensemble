function Xnew = binarizeDropMissing(X, miss)
% binarizeWithMissing  Make 1-of-K encoding of survey data.
%
%   B = binarizeDropMissing(X, miss)
%   ignores 0, which are missing values.

    [M, D] = size(X);

    dim = 1;
    for k = 1:D
        levels = unique(X(:,k))';
        for lv = levels
            if lv ~= miss
                Xnew(X(:,k) == lv, dim) = 1;
                dim = dim + 1;
            end
        end
    end

end

