function test_inference

%randn('state',5);
%rand ('state',5);

% first, make a randomly connected graph
model.nnodes   = 20;
model.ncliques = 40;
model.nvals    = 4;
rho            = 1;
maxiter        = 100;

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
psi_ij = randn(model.nvals^2,model.ncliques);
psi_i  = randn(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(psi_i,2) psi_i(:,i) = psi_i(:,i)/sum(psi_i(:,i)); end


fprintf('Testing that trw and trw_fast give the same results...\n');
tic
[b_ij b_i A]   = trw     (model,psi_ij,psi_i,rho,maxiter);
time_trw = toc;
tic
[b_ij2 b_i2 A2] = trw_fast(model,psi_ij,psi_i,rho,maxiter);
time_trw_fast = toc;
fprintf('[Times]:\ntrw:      %f\ntrw_fast: %f  \n', time_trw, time_trw_fast);
fprintf('[Difference of marginals] (should be zero) \nunivariate: %f \nbivariate:  %f \n',...
    norm(b_i-b_i2), norm(b_ij-b_ij2))
fprintf('[Approx partition functions] \ntrw:      %f \ntrw_fast: %f \ndiff:      %f (should be zero) \n',...
    A,A2,A-A2);
fprintf('\n')

fprintf('testing that meanfield and meanfield_fast give the same results...\n');
tic
[b_ij3 b_i3 A3] = meanfield(model,psi_ij,psi_i,maxiter,0);
time_meanfield = toc;
tic
[b_ij4 b_i4 A4] = meanfield_fast(model,psi_ij,psi_i,maxiter,0);
time_meanfield_fast = toc;
fprintf('[Times]\nmeanfield:      %f\nmeanfield_fast: %f  \n', time_meanfield, time_meanfield_fast);
fprintf('[Difference of marginals] (should be zero) \nunivariate: %f \nbivariate:  %f \n',...
    norm(b_i3-b_i4), norm(b_ij3-b_ij4))
fprintf('[Approx partition functions] \nmeanfield:      %f \nmeanfield_fast: %f \ndiff:             %f (should be zero) \n',...
    A3,A4,A3-A4);
fprintf('\n')

%figure(1)
subplot(2,1,1)
plots(b_i(:),b_i3(:))
title('univariate marginals')
legend('trw','meanfield')

subplot(2,1,2)
plots(b_ij(:),b_ij3(:))
title('bivariate marginals')
legend('trw','meanfield')

set(gcf,'Position',[1         400        1200         400])

return

eval_loss = @eval_ulike;

% test em likelihood thing
x = ceil(model.nvals*rand(model.nnodes/2,1));
[L dpsi_ij dpsi_i] = eval_loss(model, psi_ij, psi_i, rho, maxiter, x);
e = 1e-8;
for i=1:size(psi_i,1)
    for j=1:size(psi_i,2)
        psi_i2       = psi_i;
        psi_i2(i,j)  = psi_i(i,j) + e;
        L2 = eval_loss(model, psi_ij, psi_i2, rho, maxiter, x);
        dpsi_i2(i,j) = (1/e)*(L2-L);    
    end
end

%dpsi_i
%dpsi_i2

for i=1:size(psi_ij,1)
    for j=1:size(psi_ij,2)
        psi_ij2       = psi_ij;
        psi_ij2(i,j)  = psi_ij(i,j) + e;
        L2 = eval_loss(model, psi_ij2, psi_i, rho, maxiter, x);
        dpsi_ij2(i,j) = (1/e)*(L2-L);
    end
end

%dpsi_ij
%dpsi_ij2

fprintf('differences: %f | %f \n', norm(dpsi_i-dpsi_i2), norm(dpsi_ij-dpsi_ij2));
keyboard


v_i  = randn(size(b_i ));
v_ij = 0*randn(size(b_ij));
function [L db_ij db_i] = myloss(b_ij, b_i)
    L    = sum(v_i(:).*b_i(:)) + sum(v_ij(:).*b_ij(:));
    db_i  = v_i;
    db_ij = v_ij;
end

%e = 1e-5;
% [b_ij2 b_i2] = trw_fast(model,exp(log(psi_ij)+e*db_ij), exp(log(psi_i)+e*db_i), rho, maxiter);
% dlogpsi_ij = (1/e)*(b_ij2-b_ij);
% dlogpsi_i  = (1/e)*(b_i2 -b_i );
% 
% dpsi_ij = dlogpsi_ij./psi_ij;
% dpsi_i  = dlogpsi_i ./psi_i ;

[L dpsi_ij dpsi_i] = implicit_diff(model, psi_ij, psi_i, rho, maxiter, @myloss)

e = 1e-6;
for i=1:size(psi_i,1)
    for j=1:size(psi_i,2)
        psi_i2 = psi_i;
        psi_i2(i,j) = psi_i(i,j) + e;
        [b_ij b_i] = trw_fast(model,psi_ij,psi_i2,rho,maxiter);
        L2    = myloss(b_ij, b_i);
        dpsi_i2(i,j) = (1/e)*(L2-L);    
    end
end

%dpsi_i
%dpsi_i2

for i=1:size(psi_ij,1)
    for j=1:size(psi_ij,2)
        psi_ij2 = psi_ij;
        psi_ij2(i,j) = psi_ij(i,j) + e;
        [b_ij b_i] = trw_fast(model,psi_ij2,psi_i,rho,maxiter);
        L2    = myloss(b_ij, b_i);
        
        dpsi_ij2(i,j) = (1/e)*(L2-L);
    end
end

%dpsi_ij 
%dpsi_ij2

max(abs(dpsi_i(:)  - dpsi_i2(:)))
max(abs(dpsi_ij(:) - dpsi_ij2(:)))

figure(1)
imagesc([dpsi_i dpsi_i2]); colorbar
figure(2)
imagesc([dpsi_ij dpsi_ij2]); colorbar
end