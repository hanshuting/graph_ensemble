function [L b_ij b_i dF dG] = crf_linear_linear(model,F,G,x,z,y,rho,loss_spec)

% calculates the loss and derivative for a CRF of the form
% p(y|x) = exp( \sum_i  x(i,:) * F(:,y(i))
%               + \sum_ij z(ij,:) * G(ij,y(i)y(j))
%               - A(x))

% description of inputs:
%   model is a model
%   y is a nnodes x 1 vector of integers
%   x is a nnodes x nvals      matrix of doubles
%   z is a npairs x nvals^2    matrix of doubles
%
%   *KT* correction:
%   x is a nnodes x nfeat      matrix of doubles
%   z is a npairs x nefeat     matrix of doubles
%   G is a nvals*nvals x nefeat matrix of reals (don't have to be positive)
%
%
%   F is a nvals  x nfeat       matrix of reals  (don't have to be positive)
%   G is a nfeat  x nvals*nvals matrix of reals  (don't have to be positive)
%   rho is a double (edge appearance probability)
%   lossname is the name of a loss ('ul','cl','em', whatever)

% this is 'linear' since the univariate potential depends linearly on x
% this is 'linear' since the bivariate potentials are linear on z

% first, create psi_i and psi_ij
theta_i  = F*x';
%psi_ij = repmat( exp(G(:)), 1, model.ncliques);
theta_ij = G*z';

% [L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
%     rho, maxiter, convthresh, y, lossname, dorec);
if nargout < 4
    [L b_ij b_i] = loss_dispatch(model, theta_ij, theta_i,...
        rho, loss_spec, y);
    return
else
[L b_ij b_i dtheta_ij dtheta_i] = loss_dispatch(model, theta_ij, theta_i,...
    rho, loss_spec, y);
end

% now, need to push back derivs to input parameters

% these come from le vector chain rule
dF = dtheta_i*x;
%dG = reshape(sum(dpsi_ij.*psi_ij,2),model.nvals,model.nvals);
dG = dtheta_ij*z;
