function [L b_ij b_i dtheta_ij dtheta_i] = em(model, theta_ij, theta_i, rho, ...
    maxiter, convthresh, x, damp, inference)

% runs em with inference run to completion (non-truncated)

if nargin < 8
    damp = 0;
end

inference_stuff = regexp(inference,'/','split');
if length(inference_stuff)==1
    inference1 = inference;
    inference2 = inference;
else
    inference1 = inference_stuff{1};
    inference2 = inference_stuff{2};
end

% x is a vector of observed values
%
% zeros correspond to hidden variables
% 
% if length(x) < model.nnodes this function evaluates the 
% "approximate variational generalized em" (ha!)
% with the first length(x) variables observed, the rest considered hidden.

% initial beliefs
if strcmp(inference1,'trw')
    [b_ij b_i A] = trw_fast(model,theta_ij,theta_i,rho,maxiter,damp,convthresh);
elseif strcmp(inference1,'mnf')
    [b_ij b_i A] = meanfield_fast(model,theta_ij,theta_i,maxiter,convthresh);
else
    error('unsupported inference method: %s', inference);
end

% now, clamp
goods = find(x(:));
%psi_i0 = psi_i;
%psi_i(:,goods)=0+1e-300; % this is not a perfect solution...
%psi_i(x(goods)+(goods-1)*model.nvals)=1;
theta_i0 = theta_i;
theta_i(:,goods)=-700; % exp(-700) pushes IEEE to the brink
theta_i(x(goods)+(goods-1)*model.nvals)=0;

if strcmp(inference2,'trw')
    [b2_ij b2_i A2] = trw_fast(model,theta_ij,theta_i,rho,maxiter,damp,convthresh);
    % A2 from inference is incorrect due to changed psi_i
    who_i  = b2_i >0;
    who_ij = b2_ij>0;
    noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
    noccur = repmat(noccur',model.nvals,1);
    A2 = sum(theta_ij(who_ij) .* b2_ij(who_ij)) + ...
         sum(theta_i0(who_i ) .* b2_i (who_i )) - ...
         sum((1-rho*noccur(who_i)).*b2_i(who_i).*log(b2_i(who_i))) - ...
         rho*sum(b2_ij(who_ij).*log(b2_ij(who_ij)));
    % *KT* check the size of "crap", e.g. nonzeroness of entropy at the
    % clamped values (due to the much larger parameter values arising from
    % changed regularization scale.)
   crap =          sum((1-rho*noccur(who_i)).*b2_i(who_i).*log(b2_i(who_i))) - ...
        rho*sum(b2_ij(who_ij).*log(b2_ij(who_ij)));
    if abs(crap) > 1e-12
        warning('Crap was big: %g', crap);
    end
elseif strcmp(inference2,'mnf')
    [b2_ij b2_i A2] = meanfield_fast(model,theta_ij,theta_i,maxiter,convthresh);
    
    % A2 from inference is incorrect due to changed psi_i
    who_i  = b2_i >0;
    who_ij = b2_ij>0;
    A2 = sum(b2_ij(who_ij).*theta_ij(who_ij)) + ...
         sum(b2_i (who_i ).*theta_i0(who_i))  - ...
         sum(b2_i (who_i ).*log(  b2_i(who_i)));
    
     % same but slower!
%       A2 = 0;
%       for i=1:model.nnodes
%           for xi=1:model.nvals
%               A2 = A2 + b2_i(xi,i)*log(psi_i0(xi,i));
%           end
%           if ~x(i)
%               A2 = A2 - b2_i(x(i),i)*log(  b2_i(x(i),i));
%           end
%       end
%       for n=1:model.ncliques
%           i = model.pairs(n,1);
%           j = model.pairs(n,2);
%           who = x(i) + (x(j)-1)*model.nvals;
%           %if ~x(i) || ~x(j)
%           %    'ghgjd'
%           for who=1:model.nvals^2
%               A2 = A2 + b2_ij(who,n)*log(psi_ij(who,n));
%           end
%       end
else
    error('unsupported inference method: %s', inference);
end

L          = -A2    + A;
dtheta_ij = -b2_ij + b_ij;
dtheta_i  = -b2_i  + b_i;
%dpsi_ij    = dlogpsi_ij ./ psi_ij;
%dpsi_i     = dlogpsi_i  ./ psi_i0;
