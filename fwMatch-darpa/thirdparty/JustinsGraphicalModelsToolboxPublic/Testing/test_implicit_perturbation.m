function test_trw_bprop

randn('state',5);
rand ('state',5);

% first, make a randomly connected graph
model.nnodes   = 20;
model.ncliques = 40;
model.nvals    = 2;
rho            = .5;
maxiter        = 25;

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
    where=where+1;
end

%model.pairs

% % compute neighborhoods
% model.N1 = cell(model.nnodes,1);
% model.N2 = cell(model.nnodes,1);
% for i=1:model.nnodes
%     model.N1{i} = [];
%     model.N2{i} = [];
% end
% for c=1:model.ncliques
%     i = model.pairs(c,1);
%     j = model.pairs(c,2);
%     model.N1{i} = [model.N1{i} c];
%     model.N2{j} = [model.N2{j} c];
% end

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
psi_ij = rand(model.nvals^2,model.ncliques);
psi_i  = rand(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(psi_i,2) psi_i(:,i) = psi_i(:,i)/sum(psi_i(:,i)); end

    % you want a random loss function!?  Here's a random loss function.
    function [L db_i db_ij] = loss(b_i, b_ij)
        L     = sum(b_i(:).^2) + sum(b_ij(:).^2);
        db_i  = 2*b_i;
        db_ij = 2*b_ij;
    end

fprintf('running TRW backprop with very tight convergence threshold iterations (fast)...\n')
[L  b_ij  b_i  dpsi_ij  dpsi_i ]   = trw_bprop_fast(model,psi_ij,psi_i,rho,1e5,1e-7,@loss,0);
fprintf('running implicit with a convthresh of 1e-10...\n')
[L2 b_ij2 b_i2 dpsi_ij2 dpsi_i2]   = implicit_perturbation(model,psi_ij,psi_i,rho,100,@loss,1e-10,2,'trw');

fprintf('univariate differences of bprop_fast and implicit gradient (should be zero): %f \n', ...
mean(abs(dpsi_i(:)-dpsi_i2(:))));
fprintf('bivariate  differences of bprop_fast and implicit gradient (should be zero): %f \n', ...
mean(abs(dpsi_ij(:)-dpsi_ij2(:))));

fprintf('\n');

fprintf('running meanfield backprop with very tight convergence threshold iterations (fast)...\n')
[L  b_ij  b_i  dpsi_ij  dpsi_i ]   = meanfield_bprop_fast(model,psi_ij,psi_i,1e5,1e-7,@loss,0);
fprintf('running implicit with a convthresh of 1e-10...\n')
[L2 b_ij2 b_i2 dpsi_ij2 dpsi_i2]   = implicit_perturbation(model,psi_ij,psi_i,[],100,@loss,1e-10,2,'mnf');

fprintf('univariate differences of bprop_fast and implicit gradient (should be zero): %f \n', ...
mean(abs(dpsi_i(:)-dpsi_i2(:))));
fprintf('bivariate  differences of bprop_fast and implicit gradient (should be zero): %f \n', ...
mean(abs(dpsi_ij(:)-dpsi_ij2(:))));


end