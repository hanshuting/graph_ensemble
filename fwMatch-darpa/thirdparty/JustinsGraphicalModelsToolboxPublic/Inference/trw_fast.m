function [b_ij b_i A Z_ij reps] = trw_fast(model,theta_ij,theta_i,rho,maxiter,damp,convthresh)

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

if nargin < 6 || isempty(damp)
    damp = 0;
end
if nargin < 7
    %convthresh = .000002; % ~ 0.001 accuracy
    convthresh = .00002; % ~ 0.01 accuracy
end
%convthresh = 0;

% initialize messages
n  = ones(model.nvals,model.nnodes  );
m1 = ones(model.nvals,model.ncliques); % messages to first variable in clique
m2 = ones(model.nvals,model.ncliques); % messages to second variable

psi_ij = exp(theta_ij);
psi_i  = exp(theta_i);
psi_ij_rho = psi_ij.^(1/rho);

b_i  = n + 1e-100; % force deep copy
b_ij = zeros(size(psi_ij));
Z_ij = zeros(1,model.ncliques);

reps = trw_helper(model, psi_ij_rho, psi_i, rho, maxiter, damp, convthresh, n, m1, m2, b_i, b_ij, Z_ij);

who_i  = b_i >0;
who_ij = b_ij>0;

% count number of times each variable occurs in a clique
noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
noccur = repmat(noccur',model.nvals,1);

%who2_ij = (psi_ij>0);
%who2_i  = (psi_i >0);

A = sum(theta_ij(who_ij).*b_ij(who_ij)) + ...
    sum(theta_i (who_i ).*b_i (who_i )) - ...
    sum((1-rho*noccur(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
    rho*sum(b_ij(who_ij).*log(b_ij(who_ij)));