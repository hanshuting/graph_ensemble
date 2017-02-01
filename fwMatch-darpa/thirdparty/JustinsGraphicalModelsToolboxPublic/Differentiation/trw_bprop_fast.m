%function [L b_ij b_i dpsi_ij dpsi_i Z_ij] = trw_bprop_fast(model,psi_ij,psi_i,rho,maxiter,loss,mode)
function [L b_ij b_i dtheta_ij dtheta_i] = trw_bprop_fast(model,theta_ij,theta_i,rho,...
    maxiter,convthresh,loss,dorec)

% do truncated or non-truncated back-belief propagation

% loss is a function object that will be called like so:
% [L db_i db_ij] = loss(b_i, b_ij);

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

debug = 0;

if nargin < 8
    dorec = 1; % means do record messages
end

% initialize messages
n  = ones(model.nvals,model.nnodes  );
m1 = ones(model.nvals,model.ncliques); % messages to first variable in clique
m2 = ones(model.nvals,model.ncliques); % messages to second variable

% all messages are stored here before being updated
if dorec
    wpred = 2*maxiter*model.ncliques*model.nvals;
    mstor = zeros(wpred,1); w = 0;
else
    w     = [];
    mstor = [];
end

b_i = n;

psi_i   = exp(theta_i);
psi_ij0 = exp(theta_ij);

psi_ij = psi_ij0.^(1/rho);
b_ij0  = psi_ij; b_ij0(1)=b_ij0(1)+1e-100;

if min(min(psi_i(:)),min(psi_ij0(:)))==0 || max(max(psi_i(:)),max(psi_ij0(:)))==inf
    L = 1e9;
    b_ij = b_ij0;
    dtheta_ij = zeros(size(theta_ij));
    dtheta_i  = zeros(size(theta_i));
    return
end

if debug, tic; end
%trw_helper(model, psi_ij_rho, psi_i, rho, maxiter, damp, convthresh, n, m1, m2, b_i, b_ij, Z_ij);
actualiters = 0;
w           = 0;

trw_bprop_helper1(model, psi_ij, psi_i, rho, maxiter, convthresh, ...
    n, m1, m2, b_i, mstor, b_ij0, dorec, actualiters, w);
maxiter = actualiters;

% convthresh is not looked at again
%clear convthresh
mstor = mstor(1:w); % only use messages actually filled

if debug, time_fwprop = toc, end
if debug, tic; end

% for c=1:model.ncliques
%     i = model.pairs(c,1);
%     j = model.pairs(c,2);
%     for yi=1:model.nvals
%         for yj=1:model.nvals
%             index = yi + (yj-1)*model.nvals;
%             b_ij0(index,c) = b_ij0(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
%         end
%     end
% end

% c=1:model.ncliques;
% i = model.pairs(c,1);
% j = model.pairs(c,2);
% for yi=1:model.nvals
%     for yj=1:model.nvals
%         index = yi + (yj-1)*model.nvals;
%         b_ij0(index,c) = b_ij0(index,c).*psi_i(yi,i).*psi_i(yj,j).*n(yi,i).*n(yj,j)./m1(yi,c)./m2(yj,c);
%     end
% end

Z_ij = sum(b_ij0,1);
Z_ij = repmat(Z_ij,model.nvals^2,1);
b_ij = b_ij0 ./ Z_ij;

b_i0 = n.*psi_i;
Z_i  = sum(b_i0,1);
Z_i  = repmat(Z_i,model.nvals,1);
b_i  = b_i0 ./ Z_i;

if debug, time_int1 = toc, end

if debug, tic; end
% compute loss
[L db_i db_ij] = loss(b_i, b_ij);
if nargout<=3, return; end
% propagate back to messages
dm1  = 0*m1;
dm2  = 0*m2;

% get derivs w.r.t. unnormalized beliefs
db_i0  = db_i  ./ Z_i  - repmat(sum(db_i .*b_i0 ./Z_i .^2,1),model.nvals  ,1);
db_ij0 = db_ij ./ Z_ij - repmat(sum(db_ij.*b_ij0./Z_ij.^2,1),model.nvals^2,1);

% starter propagations (onto params and n)
dn      = db_i0 .* psi_i;
dpsi_i  = db_i0 .*n;
dpsi_ij = db_ij0.*b_ij0./psi_ij;

if debug, time_int2 = toc, end

% if debug, tic; end
% for c=1:model.ncliques
%     i = model.pairs(c,1);
%     j = model.pairs(c,2);
%     for yi=1:model.nvals
%         for yj=1:model.nvals
%             index = yi + (yj-1)*model.nvals;
%             %b_ij0(index,c) = b_ij0(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
%             dpsi_i(yi,i) = dpsi_i(yi,i) + db_ij0(index,c)*b_ij0(index,c)/psi_i(yi,i);
%             dpsi_i(yj,j) = dpsi_i(yj,j) + db_ij0(index,c)*b_ij0(index,c)/psi_i(yj,j);
%             dn(yi,i)     = dn(yi,i)     + db_ij0(index,c)*b_ij0(index,c)/n(yi,i);
%             dn(yj,j)     = dn(yj,j)     + db_ij0(index,c)*b_ij0(index,c)/n(yj,j);
%             dm1(yi,c)    = dm1(yi,c)    - db_ij0(index,c)*b_ij0(index,c)/m1(yi,c);
%             dm2(yj,c)    = dm2(yj,c)    - db_ij0(index,c)*b_ij0(index,c)/m2(yj,c);
%         end
%     end
% end
% if debug, time_int3 = toc, end

% tic
% trw_bprop_helper2(model, psi_ij, psi_i, rho, maxiter, n, m1, m2, b_i, mstor,...
%                   dm1, dm2, dn, dpsi_i, dpsi_ij, b_ij0, db_ij0, dorec);
% toc
% keyboard

if debug, tic; end
trw_bprop_helper2(model, psi_ij, psi_i, rho, maxiter, n, m1, m2, b_i, mstor,...
                  dm1, dm2, dn, dpsi_i, dpsi_ij, b_ij0, db_ij0, dorec);
if debug, time_brop = toc, end

% psi_ij = psi_ij.^(1/rho);
dpsi_ij = (1/rho)*psi_ij.^(1-rho).*dpsi_ij;

dtheta_i  = dpsi_i .*psi_i;
dtheta_ij = dpsi_ij.*psi_ij0;

% if isnan(sum(dpsi_i(:))+sum(dpsi_ij(:)))
%     keyboard
% end