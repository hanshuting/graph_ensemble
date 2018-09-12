function [L b_ij b_i dtheta_ij dtheta_i reps] = implicit_perturbation(model,theta_ij,theta_i,rho,maxiter,loss,convthresh,mode,inference,esize)

% implement the main result from Implicit Differentiation by Perturbation
% this is all a little wierd because psi_ij and psi_i are parameters
% AFTER exponentiation, whereas the paper calculates derivs before.

% loss is a function object that will be called like so:
% [L db_i db_ij] = loss(b_i, b_ij);

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

if nargin < 8
    % two sided differences are the default
    mode = 2;
end

% initial beliefs
if strcmp(inference,'trw')
    [b_ij b_i ~, ~, reps] = trw_fast(model,theta_ij,theta_i,rho,maxiter,0,convthresh);
elseif strcmp(inference,'mnf')
    [b_ij b_i A reps] = meanfield_fast(model,theta_ij,theta_i,maxiter, convthresh);
else
    error('only support mnf or trw inference');
end

if reps==maxiter
    warning('in implicit_perturbation, reps reached maxiter.');
end

% calculate loss
[L db_i db_ij] = loss(b_i, b_ij);

% calculate step size
l_inf_theta = max(max(abs(theta_ij(:))), max(abs(theta_i(:))));
l_inf_db    = max(max(abs(   db_ij(:))), max(abs(   db_i(:))));
e = sqrt(eps)*(1+l_inf_theta)/l_inf_db;

if ~exist('esize','var')
    esize = 1;
end
e = e*esize;

if strcmp(inference,'trw')
    if mode==1
        [b_ij_pos b_i_pos ~, ~, reps_pos] = trw_fast(model,theta_ij+e*db_ij,theta_i+e*db_i,rho,maxiter,0,convthresh);
        dtheta_ij = 1/(e)*(b_ij_pos - b_ij);
        dtheta_i  = 1/(e)*(b_i_pos  - b_i);
        reps = reps + reps_pos;
    elseif mode==2
        [b_ij_pos b_i_pos ~, ~, reps_pos] = trw_fast(model,theta_ij+e*db_ij,theta_i+e*db_i,rho,maxiter,0,convthresh);
        [b_ij_neg b_i_neg ~, ~, reps_neg] = trw_fast(model,theta_ij-e*db_ij,theta_i-e*db_i,rho,maxiter,0,convthresh);
        dtheta_ij = 1/(2*e)*(b_ij_pos - b_ij_neg);
        dtheta_i  = 1/(2*e)*(b_i_pos  - b_i_neg );
        reps = reps + reps_pos + reps_neg;
        
    elseif mode==4
        [b_ij_pos1 b_i_pos1 ~, ~, reps_pos1] = trw_fast(model,theta_ij+  e*db_ij,theta_i+  e*db_i,rho,maxiter,0,convthresh);
        [b_ij_pos2 b_i_pos2 ~, ~, reps_pos2] = trw_fast(model,theta_ij+2*e*db_ij,theta_i+2*e*db_i,rho,maxiter,0,convthresh);
        [b_ij_neg1 b_i_neg1 ~, ~, reps_neg1] = trw_fast(model,theta_ij-  e*db_ij,theta_i-  e*db_i,rho,maxiter,0,convthresh);
        [b_ij_neg2 b_i_neg2 ~, ~, reps_neg2] = trw_fast(model,theta_ij-2*e*db_ij,theta_i-2*e*db_i,rho,maxiter,0,convthresh);
        dtheta_ij = 1/(12*e)*(-b_ij_pos2 + 8*b_ij_pos1 - 8*b_ij_neg1 + b_ij_neg2);
        dtheta_i  = 1/(12*e)*(-b_i_pos2  + 8*b_i_pos1  - 8*b_i_neg1  + b_i_neg2 );
        reps = reps + reps_pos1 + reps_pos2 + reps_neg1 + reps_neg2;
    else
        error('currently support 1, 2, and 4 sided differences  (given %d)',mode);
    end
elseif strcmp(inference,'mnf')
    if mode==1
        [b_ij_pos b_i_pos] = meanfield_fast(model,theta_ij+e*db_ij,theta_i+e*db_i,maxiter,convthresh);
        dtheta_ij = 1/(e)*(b_ij_pos - b_ij);
        dtheta_i  = 1/(e)*(b_i_pos  - b_i);
    elseif mode==2
        [b_ij_pos b_i_pos] = meanfield_fast(model,theta_ij+e*db_ij,theta_i+e*db_i,maxiter,convthresh);
        [b_ij_neg b_i_neg] = meanfield_fast(model,theta_ij-e*db_ij,theta_i-e*db_i,maxiter,convthresh);
        dtheta_ij = 1/(2*e)*(b_ij_pos - b_ij_neg);
        dtheta_i  = 1/(2*e)*(b_i_pos  - b_i_neg );
    elseif mode==4
        [b_ij_pos1 b_i_pos1] = meanfield_fast(model,theta_ij+  e*db_ij,theta_i+  e*db_i,maxiter,convthresh);
        [b_ij_pos2 b_i_pos2] = meanfield_fast(model,theta_ij+2*e*db_ij,theta_i+2*e*db_i,maxiter,convthresh);
        [b_ij_neg1 b_i_neg1] = meanfield_fast(model,theta_ij-  e*db_ij,theta_i-  e*db_i,maxiter,convthresh);
        [b_ij_neg2 b_i_neg2] = meanfield_fast(model,theta_ij-2*e*db_ij,theta_i-2*e*db_i,maxiter,convthresh);
        dtheta_ij = 1/(12*e)*(-b_ij_pos2 + 8*b_ij_pos1 - 8*b_ij_neg1 + b_ij_neg2);
        dtheta_i  = 1/(12*e)*(-b_i_pos2  + 8*b_i_pos1  - 8*b_i_neg1  + b_i_neg2 );
    else
        error('currently support 1, 2, and 4 sided differences');
    end
    reps = 0;
else
    error('only support mnf or trw inference');
end