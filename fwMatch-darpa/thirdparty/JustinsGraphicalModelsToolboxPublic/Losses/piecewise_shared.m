function [L b_ij b_i dtheta_ij dtheta_i] = piecewise_shared(model, theta_ij, theta_i, x)

b_ij    = nan + 0*theta_ij;
b_i     = nan + 0*theta_i;

L          = 0;
dtheta_ij = 0*theta_ij;
dtheta_i  = 0*theta_i;

for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    if x(i)==0 || x(j)==0, continue; end
    
    r = .5;
    
    % compute probability distribution over single pair
    phip = zeros(model.nvals,model.nvals);
    for xi=1:model.nvals
        for xj=1:model.nvals
            who = xi + model.nvals*(xj-1);
            phip(xi,xj) = theta_ij(who,c) + r*theta_i(xi,i) + r*theta_i(xj,j);
        end
    end
    p_ij = exp(phip - log_sum_exp(phip(:)));
    p_i  = sum(p_ij,2);
    p_j  = sum(p_ij,1);
    phip = phip(:);
    p_ij = p_ij(:);
    
    who = x(i) + model.nvals*(x(j)-1);
    L   = L - phip(who) + log_sum_exp(phip);
    
    for xi=1:model.nvals
        for xj=1:model.nvals
            who = xi + model.nvals*(xj-1);
            if xi==x(i) && xj==x(j)
                dtheta_ij(who,c) = dtheta_ij(who,c)-1;
            end
            dtheta_ij(who,c) = dtheta_ij(who,c)+p_ij(who);
        end
    end
    for xi=1:model.nvals
        if xi==x(i)
            dtheta_i(xi,i) = dtheta_i(xi,i)-r*1;
        end
        dtheta_i(xi,i) = dtheta_i(xi,i) + r*p_i(xi);
    end
    for xj=1:model.nvals
        if xj==x(j)
            dtheta_i(xj,j) = dtheta_i(xj,j)-r*1;
        end
        dtheta_i(xj,j) = dtheta_i(xj,j) + r*p_j(xj);
    end
end
