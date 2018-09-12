function [L b_ij b_i dtheta_ij dtheta_i] = pseudo(model, theta_ij, theta_i, x)

b_ij    = nan + 0*theta_ij;
b_i     = nan + 0*theta_i;

L       = 0;
dtheta_ij = 0*theta_ij;
dtheta_i  = 0*theta_i;

for i=1:model.nnodes
    if x(i)==0, continue; end
    % get conditional distribution
    logb = zeros(model.nvals,1);
    bad  = 0; % flag for hidden variables
    for xi=1:model.nvals
        logb(xi) = theta_i(xi,i);
        
        for c=model.N1(i,:)
            if c==-1, break; end
            j = model.pairs(c,2);
            xj = x(j);
            if xj==0, bad=1; continue; end
            where = xi + (xj-1)*model.nvals;
            
            logb(xi) = logb(xi) + theta_ij(where,c);
        end
        for c=model.N2(i,:)
            if c==-1, break; end
            j = model.pairs(c,1);
            xj = x(j);
            if xj==0, bad=1; continue; end
            where = xj + (xi-1)*model.nvals;
            
            logb(xi) = logb(xi) + theta_ij(where,c);
        end
    end
    
    if bad, continue; end
    L = L - logb(x(i)) + log_sum_exp(logb);
    
    b = exp(logb - log_sum_exp(logb));
    %dL/dt(xi) = I[xi==x(i)] - b(xi)
    for xi=1:model.nvals
        mult = double(xi==x(i)) - b(xi);
        
        dtheta_i(xi,i) = dtheta_i(xi,i) - mult;
        
        for c=model.N1(i,:)
            if c==-1, break; end
            j = model.pairs(c,2);
            xj = x(j);
            where = xi + (xj-1)*model.nvals;
            
            dtheta_ij(where,c) = dtheta_ij(where,c) - mult;
        end
        for c=model.N2(i,:)
            if c==-1, break; end
            j = model.pairs(c,1);
            xj = x(j);
            where = xj + (xi-1)*model.nvals;
            
            dtheta_ij(where,c) = dtheta_ij(where,c) - mult;
        end
    end
end