function [L b_ij b_i dF dG] = crf_linear_linear_noties(model,F,G,x,z,y,rho,loss_spec)

% calculates the loss and derivative for a CRF of the form
% p(y|x) = exp( \sum_i  x(i,:) * F(:,y(i))
%               + \sum_ij z(ij,:) * G(ij,y(i)y(j))
%               - A(x))

% description of inputs:
%   model is a model
%   y is a nnodes x 1 vector of integers
%   x is a nnodes x nvals                  matrix of doubles
%   z is a nnodes x nvals^2                matrix of doubles
%   F is a nvals x nfeat x nnodes          matrix of reals  (don't have to be positive)
%   G is a nfeat x nvals*nvals x ncliques  matrix of reals  (don't have to be positive)
%   rho is a double (edge appearance probability)
%   lossname is the name of a loss ('ul','cl','em', whatever)

% this is 'linear' since the univariate potential depends linearly on x
% this is 'independent' since the bivariate potentials are independent of x

[nvals nfeats nnodes] = size(F);
[~,~,ncliques]        = size(G);

% first, create psi_i and psi_ij
% psi_i = zeros(nvals,nnodes);
% for i=1:nnodes
%     psi_i(:,i)  = exp( F(:,:,i)*x(i,:)');
% end
% works but only faster on big matrices
theta_i = multiprod(F,x',[1 2],1);

% psi_ij = zeros(nvals^2,ncliques);
% for c=1:ncliques
%     psi_ij(:,c) = exp( G(:,:,c)*z(c,:)');
% end
theta_ij = squeeze(multiprod(G,z',[1 2],1));

[L b_ij b_i dtheta_ij dtheta_i] = loss_dispatch(model, theta_ij, theta_i,...
    rho, loss_spec, y);

% now, need to push back derivs to input parameters

% these come from le vector chain rule
% dF = 0*F;
% for i=1:nnodes
%     dF(:,:,i) = (dpsi_i(:,i).*psi_i(:,i))*x(i,:);
% end

%dF = multiprod(dpsi_i.*psi_i,reshape(x',[nfeats,1,size(x,1)]),1,2);
dF = multiprod(dtheta_i,reshape(x',[nfeats,1,size(x,1)]),1,2);
dF = permute(dF,[2 1 3]);

% dG = 0*G;
% for c=1:ncliques
%     dG(:,:,c) = (dpsi_ij(:,c).*psi_ij(:,c))*z(c,:);
% end
dG = multiprod(dtheta_ij,reshape(z',[size(z,2),1,size(z,1)]),1,2);
dG = permute(dG,[2 1 3]);