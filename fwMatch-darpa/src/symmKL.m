function d = symmKL(p, q)
% symmKL  Symmetrized KL divergence (no smoothing)
    d = sum(p .* (log(p) - log(q))) + ...
        sum(q .* (log(q) - log(p)));
    
end
