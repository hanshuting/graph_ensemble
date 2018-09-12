function [L b_ij b_i dF dG] = crf_linear_independent(model,F,G,x,y,rho,loss_spec)

% calculates the loss and derivative for a CRF of the form
% p(y|x) = exp( \sum_i  x(i,:) * F(:,y(i))
%               + \sum_ij G(y(i),y(j))
%               - A(x))

% description of inputs:
%   model is a model
%   y is a nnodes x 1     vector of integers
%   x is a nnodes x nfeat matrix of doubles
%   F is a nvals  x nfeat  matrix of reals  (don't have to be positive)
%   G is a nvals  x nvals matrix of reals  (don't have to be positive)
%   rho is a double (edge appearance probability)
%   lossname is the name of a loss ('ul','cl','em', whatever)

% this is 'linear' since the univariate potential depends linearly on x
% this is 'independent' since the bivariate potentials are independent of x

% maxiter    = 2;
% convthresh = 0;
% dorec      = 1;

% % first, create psi_i and psi_ij
% psi_i  = exp( F*x');
% psi_ij = repmat( exp(G(:)), 1, model.ncliques);
% 
% %[L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
% %    rho, maxiter, convthresh, y, lossname, dorec);
% [L b_ij b_i dpsi_ij dpsi_i] = loss_dispatch(model, psi_ij, psi_i,...
%     rho, loss_spec, y);
% 
% % now, need to push back derivs to input parameters
% 
% % these come from le vector chain rule
% dF = (dpsi_i.*psi_i)*x;
% dG = reshape(sum(dpsi_ij.*psi_ij,2),model.nvals,model.nvals);

% instead, just call the linear method

z = ones(model.ncliques,1);
[L b_ij b_i dF dG] = crf_linear_linear(model,F,G(:),x,z,y,rho,loss_spec);
dG = reshape(dG,size(G));
