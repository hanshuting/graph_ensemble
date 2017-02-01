function test_trw_bprop

randn('state',5);
rand ('state',5);

% first, make a randomly connected graph
model.nnodes   = 20;
model.ncliques = 40;
model.nvals    = 2;
rho            = .5;
%maxiter        = 500;
%convthresh     = 1e-7;
%dorec          = 0;

maxiterAR     = {500 , 10, 50};
convthreshAR  = {1e-6, 0, 1e-4};
dorecAR       = {0   , 1, 1};
descAR        = {'Back Belief Propagation', 'Truncated Fitting', ...
    'Record messages, but w/ a conv. threshold'};

% you want a random loss function!?  Here's a random loss function.
    function [L db_i db_ij] = loss(b_i, b_ij)
        L     = sum(b_i(:).^2) + sum(b_ij(:).^2);
        db_i  = 2*b_i;
        db_ij = 2*b_ij;
    end


for N=2:length(maxiterAR)
    maxiter    = maxiterAR{N};
    convthresh = convthreshAR{N};
    dorec      = dorecAR{N};
    desc       = descAR{N};
    fprintf('testing with maxiter=%d  convthresh=%f  dorec=%d  (%s)\n',...
        maxiter,convthresh,dorec,desc);
    
model = gridmodel(50,50,8);

% now, make some random potentials
psi_ij = rand(model.nvals^2,model.ncliques);
psi_i  = rand(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(psi_i,2) psi_i(:,i) = psi_i(:,i)/sum(psi_i(:,i)); end

fprintf('running backprop (fast)...\n')
tic
[L2 b_ij2 b_i2 dpsi_ij2 dpsi_i2]   = trw_bprop_scheduled_fast(model,psi_ij,psi_i,rho,maxiter,convthresh,@loss,dorec);
toc

fprintf('running backprop...\n')
tic
[L b_ij b_i dpsi_ij dpsi_i]   = trw_bprop_scheduled(model,psi_ij,psi_i,rho,maxiter,convthresh,@loss,dorec);
toc

fprintf('univariate differences of bprop and bprop_fast gradient (should be zero): %f \n', ...
mean(abs(dpsi_i(:)-dpsi_i2(:))));
fprintf('bivariate  differences of bprop and bprop_fast gradient (should be zero): %f \n', ...
mean(abs(dpsi_ij(:)-dpsi_ij2(:))));

fprintf('computing univariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-8;
for i=1:size(psi_i,1)
    for j=1:size(psi_i,2)
        psi_i2       = psi_i;
        psi_i2(i,j)  = psi_i(i,j) + e;
        L2 = trw_bprop_scheduled(model,psi_ij,psi_i2,rho,maxiter,convthresh,@loss,dorec);
        dpsi_i2(i,j) = (1/e)*(L2-L);
        where = where + 1;
        printstatus(where / numel(psi_i), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')

subplot(2,1,1); plots(dpsi_i(:),dpsi_i2(:));
title('dL/dtheta for univariate theta');

fprintf('differences of bprop and finite diff gradient (should be zero): %f \n', ...
mean(abs(dpsi_i(:)-dpsi_i2(:))));

fprintf('computing bivariate diffs by finite differences (this is slow!)...\n')
where = 1;
for i=1:size(psi_ij,1)
    for j=1:size(psi_ij,2)
        psi_ij2       = psi_ij;
        psi_ij2(i,j)  = psi_ij(i,j) + e;
        L2 = trw_bprop_scheduled(model,psi_ij2,psi_i,rho,maxiter,convthresh,@loss,dorec);
        dpsi_ij2(i,j) = (1/e)*(L2-L);
        where = where + 1;
        printstatus(where / numel(psi_ij), 50, i~=1 || j~=1 );
    end
end
fprintf('\n')
fprintf('differences of bprop and finite diff gradient (should be zero): %f \n', ...
mean(abs(dpsi_ij(:)-dpsi_ij2(:))));


subplot(2,1,2); plots(dpsi_ij(:),dpsi_ij2(:));
title('dL/dtheta for bivariate theta');

fprintf('\n\n\n\n');

end


end