function test_em_trunc

randn('state',5);
rand ('state',5);

% first, make a randomly connected graph
model.nnodes   =   10;
model.ncliques =   20;
model.nvals    =    3;
rho            =   .5;
maxiter        =    5;
convthresh     =    0;

%inference = 'trw/trw';
%inference = 'trw/mnf';
%inference = 'mnf/trw';
inference = 'mnf/mnf';

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
psi_ij = .5+rand(model.nvals^2,model.ncliques);
psi_i  = .5+rand(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(psi_i,2) psi_i(:,i) = psi_i(:,i)/sum(psi_i(:,i)); end

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
%[L b_ij b_i dpsi_ij dpsi_i]   = trw_bprop(model,psi_ij,psi_i,rho,maxiter,@loss);
[L b_ij b_i dpsi_ij dpsi_i] = em_trunc(model, psi_ij, psi_i, rho, maxiter, convthresh, x, inference);

fprintf('computing bivariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-8;
for i=1:size(psi_ij,1)
    for j=1:size(psi_ij,2)
        psi_ij2       = psi_ij;
        psi_ij2(i,j)  = psi_ij(i,j) + e;
        %L2 = trw_bprop(model,psi_ij2,psi_i,rho,maxiter,@loss);
        L2 = em_trunc(model, psi_ij2, psi_i, rho, maxiter, convthresh, x, inference);
        psi_ij3       = psi_ij;
        psi_ij3(i,j)  = psi_ij(i,j) - e;
        L3 = em_trunc(model, psi_ij3, psi_i, rho, maxiter, convthresh, x, inference);
        %dpsi_ij2(i,j) = (1/e)*(L2-L);
        dpsi_ij2(i,j) = (1/e/2)*(L2-L3);
        where = where + 1;
        printstatus(where / numel(psi_ij), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')
fprintf('differences of em and finite diff gradient (should be zero): %f \n', ...
mean(abs(dpsi_ij(:)-dpsi_ij2(:))));

subplot(2,1,2); plots(dpsi_ij(:),dpsi_ij2(:));
title('dL/dtheta for bivariate theta');

fprintf('computing univariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-5;
for i=1:size(psi_i,1)
    for j=1:size(psi_i,2)
        psi_i2       = psi_i;
        psi_i2(i,j)  = psi_i(i,j) + e;
        L2 = em_trunc(model, psi_ij, psi_i2, rho, maxiter, convthresh, x, inference);
        
        psi_i3       = psi_i;
        psi_i3(i,j)  = psi_i(i,j) - e;
        L3 = em_trunc(model, psi_ij, psi_i3, rho, maxiter, convthresh, x, inference);
        dpsi_i2(i,j) = (1/e/2)*(L2-L3);
        where = where + 1;
        printstatus(where / numel(psi_i), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')

subplot(2,1,1); plots(dpsi_i(:),dpsi_i2(:));
title('dL/dtheta for univariate theta');

fprintf('differences of em and finite diff gradient (should be zero): %f \n', ...
mean(abs(dpsi_i(:)-dpsi_i2(:))));

end