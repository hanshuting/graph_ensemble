function [L b_ij b_i dF dG] = crf_linear_independent_noties(model,F,G,x,y,rho,loss_spec)

% calculates the loss and derivative for a CRF of the form
% p(y|x) = exp( \sum_i  x(i,:) * F(:,y(i))
%               + \sum_ij G(y(i),y(j))
%               - A(x))

% instead, just call the linear method

z = ones(model.ncliques,1);
G2 = reshape(G,[model.nvals^2,1,model.ncliques]);
[L b_ij b_i dF dG] = crf_linear_linear_noties(model,F,G2,x,z,y,rho,loss_spec);
dG = reshape(dG,size(G));
