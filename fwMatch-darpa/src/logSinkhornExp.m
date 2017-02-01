function B = logSinkhornExp(A, prec)
% TODO: Additive version of sinkhorn?    
    if nargin < 2
        prec = .00001;
    end

    [B,x,y] = sinkhorn(exp(A), prec);
    
    B = log(B);

end

