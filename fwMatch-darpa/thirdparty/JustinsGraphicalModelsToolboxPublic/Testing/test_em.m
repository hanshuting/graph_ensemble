function test_em

randn('state',5);
rand ('state',5);

% first, make a randomly connected graph
model.nnodes   =   20;
model.ncliques =   40;
model.nvals    =    3;
rho            =   .5;
maxiter        = 1000;
damp           =    0;
convthresh     = 1e-4;

% 4 pairs to test
inference = 'trw/trw';
inference = 'trw/mnf';
%inference = 'mnf/trw';
%inference = 'mnf/mnf';

model.pairs = zeros(model.ncliques,2);
where=1;
for k=1:model.nnodes-1
    model.pairs(where,:) = [k k+1];
    where=where+1;
end
for where=where:model.ncliques
    k = ceil(rand*model.nnodes);
    l = k;
    while l==k
        l = ceil(rand*model.nnodes);
    end
    model.pairs(where,:) = [k l];
    %where=where+1;
end

maxN = 10;
model.N1 = -1*ones(model.nnodes,maxN);
model.N2 = -1*ones(model.nnodes,maxN);
where1 = ones(model.nnodes,1);
where2 = ones(model.nnodes,1);
for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    model.N1(i,where1(i))=c;
    model.N2(j,where2(j))=c;
    where1(i)=where1(i)+1;
    where2(j)=where2(j)+1;
end

% now, make some random potentials
theta_ij = rand(model.nvals^2,model.ncliques);
theta_i  = rand(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(theta_i,2) theta_i(:,i) = theta_i(:,i)/sum(theta_i(:,i)); end

    % you want a random loss function!?  Here's a random loss function.
    function [L db_i db_ij] = loss(b_i, b_ij)
        L     = sum(b_i(:).^2) + sum(b_ij(:).^2);
        db_i  = 2*b_i;
        db_ij = 2*b_ij;
    end

x = ceil(rand(model.nnodes,1)*model.nvals);
% make some hidden
x(ceil(rand(model.nnodes/2,1)*model.nnodes))=0;

fprintf('running em...\n')
%[L b_ij b_i dtheta_ij dtheta_i]   = trw_bprop(model,theta_ij,theta_i,rho,maxiter,@loss);
[L b_ij b_i dtheta_ij dtheta_i] = em(model, theta_ij, theta_i, rho, maxiter, convthresh, x, damp, inference);

fprintf('computing bivariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-8;
for i=1:size(theta_ij,1)
    for j=1:size(theta_ij,2)
        theta_ij2       = theta_ij;
        theta_ij2(i,j)  = theta_ij(i,j) + e;
        %L2 = trw_bprop(model,theta_ij2,theta_i,rho,maxiter,@loss);
        L2 = em(model, theta_ij2, theta_i, rho, maxiter, convthresh, x, damp, inference);
        theta_ij3       = theta_ij;
        theta_ij3(i,j)  = theta_ij(i,j) - e;
        L3 = em(model, theta_ij3, theta_i, rho, maxiter, convthresh, x, damp, inference);
        %dtheta_ij2(i,j) = (1/e)*(L2-L);
        dtheta_ij2(i,j) = (1/e/2)*(L2-L3);
        where = where + 1;
        printstatus(where / numel(theta_ij), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')
fprintf('differences of em and finite diff gradient (should be zero): %f \n', ...
mean(abs(dtheta_ij(:)-dtheta_ij2(:))));

subplot(2,1,2); plots(dtheta_ij(:),dtheta_ij2(:));
title('dL/dtheta for bivariate theta');


fprintf('computing univariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-8;
for i=1:size(theta_i,1)
    for j=1:size(theta_i,2)
        theta_i2       = theta_i;
        theta_i2(i,j)  = theta_i(i,j) + e;
        %L2 = trw_bprop(model,theta_ij,theta_i2,rho,maxiter,@loss);
        L2 = em(model, theta_ij, theta_i2, rho, maxiter, convthresh, x, damp, inference);
        dtheta_i2(i,j) = (1/e)*(L2-L);
        where = where + 1;
        printstatus(where / numel(theta_i), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')

subplot(2,1,1); plots(dtheta_i(:),dtheta_i2(:));
title('dL/dtheta for univariate theta');

fprintf('differences of em and finite diff gradient (should be zero): %f \n', ...
mean(abs(dtheta_i(:)-dtheta_i2(:))));


end