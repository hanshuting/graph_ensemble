function [L b_ij b_i dtheta_ij dtheta_i] = piecewise(model, theta_ij, theta_i, x)

b_ij    = nan + 0*theta_ij;
b_i     = nan + 0*theta_i;

L          = 0;
dtheta_ij = 0*theta_ij;
dtheta_i  = 0*theta_i;


psi_i  = exp(theta_i);
psi_ij = exp(theta_ij);

% L         = L + sum(log(sum(psi_i,1)));
% dlogpsi_i = dlogpsi_i + 

for i=1:model.nnodes
    if x(i)==0, continue; end
    
    L = L - theta_i(x(i),i);
    dtheta_i(x(i),i) = dtheta_i(x(i),i) - 1;
    
    L = L + log(sum(psi_i(:,i)));
    dtheta_i(:,i) = dtheta_i(:,i) + psi_i(:,i) / sum(psi_i(:,i));
end

for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    if x(i)==0 || x(j)==0, continue; end
    
    who = x(i) + model.nvals*(x(j)-1);
    
    L                = L               - theta_ij(who,c);
    dtheta_ij(who,c) = dtheta_ij(who,c) - 1;
    
    L = L + log(sum(psi_ij(:,c)));
    dtheta_ij(:,c) = dtheta_ij(:,c) + psi_ij(:,c) / sum(psi_ij(:,c));
end

%     
%     % get conditional distribution
%     
%     
%     logb = zeros(model.nnodes,1);
%     bad  = 0; % flag for hidden variables
%     for xi=1:model.nvals
%         logb(xi) = psi_i(xi,i);
%         
%         for c=model.N1(i,:)
%             if c==-1, break; end
%             j = model.pairs(c,2);
%             xj = x(j);
%             if xj==0, bad=1; continue; end
%             where = xi + (xj-1)*model.nvals;
%             
%             logb(xi) = logb(xi) + psi_ij(where,c);
%         end
%         for c=model.N2(i,:)
%             if c==-1, break; end
%             j = model.pairs(c,1);
%             xj = x(j);
%             if xj==0, bad=1; continue; end
%             where = xj + (xi-1)*model.nvals;
%             
%             logb(xi) = logb(xi) + psi_ij(where,c);
%         end
%     end