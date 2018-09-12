function [L b_ij b_i dF dG] = crf_linear_independent_pairtypes(model,...
    F,G,x,y,rho,loss_spec)

% linear crf with independent pairwise potentials with parameter ties
%
% calculates the loss and derivative for a CRF of the form
% p(y|x) = exp( \sum_i  x(i,:) * F(:,y(i))
%               + \sum_ij G_j(y(i),y(j))
%               - A(x))

% description of inputs:
%   model is a model
%   y is a nnodes x 1     vector of integers
%   x is a nnodes x nvals matrix of doubles
%   F is a nfeat  x nvals matrix of reals  (don't have to be positive)
%   G is a nvals*nvals x k matrix of reals  (don't have to be positive)
%   rho is a double (edge appearance probability)
%   lossname is the name of a loss ('ul','cl','em', whatever)

% uses model.pairtypes to choose which column (j) of G determines pairwise
% potential

% this is 'linear' since the univariate potential depends linearly on x
% this is 'independent' since the bivariate potentials are independent of x

%z = ones(model.ncliques,1);
npairtypes = max(model.pairtype);
% z = zeros(model.ncliques,npairtypes);
% for c=1:model.ncliques
%     z(c,model.pairtype(c))=1;
% end
subs = [(1:model.ncliques)' model.pairtype];
sz   = [model.ncliques,npairtypes];
vals = ones(model.ncliques,1);
z    = accumarray(subs,vals,sz);

[L b_ij b_i dF dG] = crf_linear_linear(model,F,G,x,z,y,rho,loss_spec);
