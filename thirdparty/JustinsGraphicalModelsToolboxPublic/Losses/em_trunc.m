function [L b_ij b_i dtheta_ij dtheta_i] = em_trunc(model, theta_ij, theta_i, rho, ...
    maxiter, convthresh, x, inference)

% runs truncated em

% x is a vector of observed values. zeros correspond to hidden variables
% should this be called "truncated approximate variational generalized em"!?

%maxiter    = 5;
dorec      = 1;
%convthresh = 0;

if convthresh ~= 0
    error('em_trunc is intended to be used with truncated fitting only (convthresh==0)')
end

noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
noccur = repmat(noccur',model.nvals,1);

inference_stuff = regexp(inference,'/','split');
if length(inference_stuff)==1
    inference1 = inference;
    inference2 = inference;
else
    inference1 = inference_stuff{1};
    inference2 = inference_stuff{2};
end

function [L db_i db_ij] = loss1(b_i, b_ij)
who_i  = b_i >0;
who_ij = b_ij>0;
L = sum(theta_ij(who_ij).*b_ij(who_ij)) + ...
    sum(theta_i (who_i ).*b_i (who_i )) - ...
    sum((1-rho*noccur(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
    rho*sum(b_ij(who_ij).*log(b_ij(who_ij)));
db_ij = 0*b_ij;
db_i  = 0*b_i;
db_ij(who_ij) = db_ij(who_ij) + theta_ij(who_ij);
db_i (who_i)  = db_i (who_i)  + theta_i (who_i);
db_i (who_i)  = db_i (who_i)  - (1+log(b_i(who_i))).*(1-rho*noccur(who_i));
db_ij(who_ij) = db_ij(who_ij) - (1+log(b_ij(who_ij)))*rho;
end

%[L1 b_ij b_i dpsi_ij1 dpsi_i1] = trw_bprop_fast(model,psi_ij,psi_i,rho,...
%    maxiter,convthresh,@loss1,dorec);
if strcmp(inference1,'trw')
    [L1 b_ij b_i dtheta_ij1 dtheta_i1] = trw_bprop_fast(model,theta_ij,theta_i,rho,...
        maxiter,convthresh,@loss1,dorec);
elseif strcmp(inference1,'trwpll')
    [L1 b_ij b_i dtheta_ij1 dtheta_i1] = trw_bprop_scheduled_fast(model,theta_ij,theta_i,rho,...
        maxiter,convthresh,@loss1,dorec);
elseif strcmp(inference1,'mnf')
    [L1 b_ij b_i dtheta_ij1 dtheta_i1] = meanfield_bprop_fast(model,theta_ij,theta_i,...
        maxiter,convthresh,@loss1,dorec);
else
    error('unsuppored inference method: %s', inference);
end

% psi_i0 = psi_i;
% who_i  = b_i >0;
% who_ij = b_ij>0;
% dpsi_ij10 = 0*psi_ij;
% dpsi_i10  = 0*psi_i;
% dpsi_ij10(who_ij) = b_ij(who_ij)./psi_ij(who_ij);
% dpsi_i10 (who_i)  = b_i (who_i) ./psi_i (who_i);
% dpsi_ij1 = dpsi_ij1 + dpsi_ij10;
% dpsi_i1  = dpsi_i1  + dpsi_i10;
% % now, clamp
% goods = find(x(:));
% psi_i(:,goods)=0+1e-300; % this is not a perfect solution...
% psi_i(x(goods)+(goods-1)*model.nvals)=1;


theta_i0 = theta_i;
who_i  = b_i >0;
who_ij = b_ij>0;
dtheta_ij10 = 0*theta_ij;
dtheta_i10  = 0*theta_i;
dtheta_ij10(who_ij) = b_ij(who_ij);
dtheta_i10 (who_i)  = b_i (who_i) ;
dtheta_ij1 = dtheta_ij1 + dtheta_ij10;
dtheta_i1  = dtheta_i1  + dtheta_i10;
% now, clamp
goods = find(x(:));
theta_i(:,goods)=-300; % this is not a perfect solution...
theta_i(x(goods)+(goods-1)*model.nvals)=0;


%     [b2_ij b2_i A2] = trw_fast(model,psi_ij,psi_i,rho,maxiter,damp,convthresh);
%     % A2 from inference is incorrect due to changed psi_i
%     who_i  = b2_i >0;
%     who_ij = b2_ij>0;
%     noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
%     noccur = repmat(noccur',model.nvals,1);
%     A2 = sum(log(psi_ij(who_ij)) .* b2_ij(who_ij)) + ...
%          sum(log(psi_i0(who_i )) .* b2_i (who_i )) - ...
%          sum((1-rho*noccur(who_i)).*b2_i(who_i).*log(b2_i(who_i))) - ...
%          rho*sum(b2_ij(who_ij).*log(b2_ij(who_ij)));

function [L db_i db_ij] = loss2(b_i, b_ij)
who_i  = logical( double(b_i >0) .* double(theta_i>-700));
who_ij = b_ij>0;
L = sum(theta_ij(who_ij).*b_ij(who_ij)) + ...
    sum(theta_i0(who_i ).*b_i (who_i )) - ...
    sum((1-rho*noccur(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
    rho*sum(b_ij(who_ij).*log(b_ij(who_ij)));
db_ij = 0*b_ij;
db_i  = 0*b_i;
db_ij(who_ij) = db_ij(who_ij) + theta_ij(who_ij);
db_i (who_i)  = db_i (who_i)  + theta_i0(who_i);
db_i (who_i)  = db_i (who_i)  - (1+log(b_i(who_i))).*(1-rho*noccur(who_i));
db_ij(who_ij) = db_ij(who_ij) - (1+log(b_ij(who_ij)))*rho;
end

% special case:
% if EVERYTHING is clamped, secondary inference can be really easy
% and there is no reason to do many iterations
if length(find(x(:))) == numel(x)
    maxiter = 1;
%     fprintf('easy inference')
% else
%     fprintf('hard inference');
end

if strcmp(inference2,'trw')
    [L2 b_ij2 b_i2 dtheta_ij2 dtheta_i2] = trw_bprop_fast(model,theta_ij,theta_i,rho,...
        maxiter,convthresh,@loss2,dorec);
elseif strcmp(inference2,'trwpll')
    [L2 b_ij2 b_i2 dtheta_ij2 dtheta_i2] = trw_bprop_scheduled_fast(model,theta_ij,theta_i,rho,...
        maxiter,convthresh,@loss2,dorec);
elseif strcmp(inference2,'mnf')
    [L2 b_ij2 b_i2 dtheta_ij2 dtheta_i2] = meanfield_bprop_fast(model,theta_ij,theta_i,...
        maxiter,convthresh,@loss2,dorec);
else
    error('unsuppored inference method: %s', inference);
end

% zero out clamped variables.
% dpsi_i2(:,goods) = 0;
% who_i  = logical( double(b_i2>0) .* double(psi_i>0));
% who_ij = b_ij2>0;
% dpsi_ij20 = 0*psi_ij;
% dpsi_i20  = 0*psi_i;
% dpsi_ij20(who_ij) = b_ij2(who_ij)./psi_ij(who_ij);
% dpsi_i20 (who_i)  = b_i2 (who_i) ./psi_i0(who_i);
% dpsi_ij2 = dpsi_ij2 + dpsi_ij20;
% dpsi_i2  = dpsi_i2  + dpsi_i20;
% L = L1-L2;
% dpsi_ij = dpsi_ij1 - dpsi_ij2;
% dpsi_i  = dpsi_i1  - dpsi_i2;

dtheta_i2(:,goods) = 0;
who_i  = logical( double(b_i2>0) .* double(theta_i>-700));
who_ij = b_ij2>0;
dtheta_ij20 = 0*theta_ij;
dtheta_i20  = 0*theta_i;
dtheta_ij20(who_ij) = b_ij2(who_ij);
dtheta_i20 (who_i)  = b_i2 (who_i) ;
dtheta_ij2 = dtheta_ij2 + dtheta_ij20;
dtheta_i2  = dtheta_i2  + dtheta_i20;
L = L1-L2;
dtheta_ij = dtheta_ij1 - dtheta_ij2;
dtheta_i  = dtheta_i1  - dtheta_i2;

end