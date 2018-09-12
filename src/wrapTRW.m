function obj = wrapTRW(y, lambda, N, L, edges, edgeRhos)
% wrapTRW  Return objective function for fmincon.
    
    rho = makeRhoOvercomplete(edgeRhos, edges, N, L);
    
    function [f, g] = o(mu)           
        assert(all(size(rho) == size(mu)));
        
        f = 0.5*lambda*sum((y - mu).^2) -     HTRW(mu, rho, L, N);
        g = -lambda*(y - mu)            - gradHTRW(mu, rho, L, N);
    end

    obj = @o;
end

