function W = learnClosedForm(mu, rho)
% learnClosedForm  Learn weight matrix W from samples X, closed form

    if nargin < 2
        % Use Bethe weights
        rho = ones(D, 1);
    end
    
    rho = rho(:);
    riPlusRj = bsxfun(@plus, rho, rho');
    
    W = log(mu) + (riPlusRj - 1) * log(1 - mu) + riPlusRj;

end

