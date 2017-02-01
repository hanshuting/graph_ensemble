function [L b_ij b_i dtheta_ij dtheta_i] = pseudo_fast(model, theta_ij, theta_i, x)

b_ij    = 0*theta_ij;
b_i     = 0*theta_i;

L       = 0;
dtheta_ij = 0*theta_ij;
dtheta_i  = 0*theta_i;

pseudo_helper(model,theta_ij,theta_i,x,L,dtheta_ij,dtheta_i);