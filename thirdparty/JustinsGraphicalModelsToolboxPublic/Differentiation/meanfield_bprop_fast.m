function [L b_ij b_i dtheta_ij dtheta_i Z_ij] = meanfield_bprop_fast(model,theta_ij,theta_i,...
    maxiter,convthresh,loss,dorec)

%[b_ij b_i A] = meanfield(model,psi_ij,psi_i,maxiter)

% you give a convergence threshold and a maximum number of iterations.
% forward propagation stops when EITHER of these is reached.
% backward propagation always uses the same number of iterations as
% forward.

% generally, one would use something like
% dorec=0, maxiter=inf, convthresh = 1e-5 -> Back belief propagation
% dorec=1, maxiter=5,   convthresh = 0    -> Truncated fitting
% but you can use both if you are a wierdo.  (In practice, one probably
% doesn't want to allow infinite iterations with BBP.)

% loss is a function object that will be called like so:
% [L db_i db_ij] = loss(b_i, b_ij);

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_ij should have size [model.nvals^2 size(model.pairs,1)]');
end

if nargin < 7
    dorec = 1; % means do record messages
end

if dorec==1 && maxiter > 50
    warning('dorec=%d with maxiter=%d',dorec,maxiter);
end

psi_i  = exp(theta_i);
psi_ij = exp(theta_ij);

% initial beliefs
b_i  = 1+0*psi_i;
b_ij = 1+0*psi_ij;

% all messages are stored here before being updated
if dorec
    wpred = maxiter*(2*model.nnodes-1)*model.nvals;
    mstor = zeros(wpred,1); w = 0;
else
    w     = [];
    mstor = [];
end

actualiters = 0;
w           = 0;
meanfield_bprop_helper1(model, psi_ij, psi_i, maxiter, convthresh, ...
    b_ij, b_i, mstor, dorec, actualiters, w);
maxiter = actualiters;

%fprintf('iters: %d\n', actualiters);

% for c=1:model.ncliques
%     for yi=1:model.nvals
%         for yj=1:model.nvals
%             i = model.pairs(c,1);
%             j = model.pairs(c,2);
%             index = yi + (yj-1)*model.nvals;
%             b_ij(index,c) = b_i(yi,i)*b_i(yj,j);
%         end
%     end
% end

[YJ YI] = meshgrid(1:model.nvals,1:model.nvals);
b_ij = b_i(YI,model.pairs(:,1)).*b_i(YJ,model.pairs(:,2));

% compute loss
[L db_i db_ij] = loss(b_i, b_ij);
% propagate back to messages

% for c=1:model.ncliques
%     for yi=1:model.nvals
%         for yj=1:model.nvals
%             i = model.pairs(c,1);
%             j = model.pairs(c,2);
%             index = yi + (yj-1)*model.nvals;
%             db_i(yi,i) = db_i(yi,i) + db_ij(index,c)*b_ij(index,c)/b_i(yi,i);
%             db_i(yj,j) = db_i(yj,j) + db_ij(index,c)*b_ij(index,c)/b_i(yj,j);
%         end
%     end
% end

% make a deep copy of beliefs before backpropagating
b_i2 = b_i; b_i2(1,1) = b_i2(1,1)+0;

dpsi_ij = 0*psi_ij;
dpsi_i  = 0*psi_i;

meanfield_bprop_helper2(model, psi_ij, psi_i, maxiter, convthresh, ...
    b_ij, b_i2, mstor, dorec, w, dpsi_ij, dpsi_i, db_ij, db_i);

dtheta_i  = dpsi_i .*psi_i;
dtheta_ij = dpsi_ij.*psi_ij;


end