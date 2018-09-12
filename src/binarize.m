function [ Xnew ] = binarize(X, discDims)
% binarize  Convert mixed data matrix X into binary data Xnew
%   Xnew = binarize(X, discDims) encodes discrete dimensions into one-hot
%   and continuous dimensions into above/below median, outputting Xnew.

    [N, Dorig] = size(X);   
    dim = 1;
    contDims = setdiff(1:Dorig, discDims);
    
    for k = contDims
        mm = median(X(:,k));
        if (sum(X(:,k)>mm-1e-4)==N)
            Xnew(:,dim) = X(:,k) < mm+1e-4;
        else
            Xnew(:,dim) = X(:,k) > mm-1e-4;
        end
        dim = dim + 1;
    end
    
    for k = discDims
        levels = unique(X(:,k))';
        for lv = levels
            Xnew(X(:,k) == lv, dim) = 1;
            dim = dim + 1;
        end
    end
    
    % Compatability
    %Xnew = double(Xnew);
end

