function test_implicit_loss

randn('state',5);
rand ('state',5);

% first, make a randomly connected graph
model.nnodes   = 50;
model.ncliques = 100;
model.nvals    = 2;
rho            = .5;
%maxiter        = 500;
%convthresh     = 1e-7;
%dorec          = 0;

maxiterAR     = {500 , 5, 50};
convthreshAR  = {1e-6, 0, 1e-4};
dorecAR       = {0   , 1, 1};
descAR        = {'Back Belief Propagation', 'Truncated Fitting', ...
    'Record messages, but w/ a conv. threshold'};

lossname = 'ul';
x        = ceil(rand(model.nnodes,1)*model.nvals);
x(2:2:end)=0;

for N=1:length(maxiterAR)
    maxiter    = maxiterAR{N};
    convthresh = convthreshAR{N};
    dorec      = dorecAR{N};
    desc       = descAR{N};
    fprintf('testing with maxiter=%d  convthresh=%f  dorec=%d  (%s)\n',...
        maxiter,convthresh,dorec,desc);

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

fprintf('running backprop...\n')
%[L b_ij b_i dpsi_ij dpsi_i]     = trw_bprop(model,psi_ij,psi_i,rho,maxiter,convthresh,@loss,dorec);
[L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model,psi_ij,psi_i,rho,maxiter,convthresh,x,lossname,dorec);


fprintf('computing univariate diffs by finite differences (this is slow!)...\n')
where = 1;
e = 1e-8;
for i=1:size(psi_i,1)
    for j=1:size(psi_i,2)
        psi_i2       = psi_i;
        psi_i2(i,j)  = psi_i(i,j) + e;
        %L2 = trw_bprop(model,psi_ij,psi_i2,rho,maxiter,convthresh,@loss,dorec);
        L2 = implicit_loss(model,psi_ij,psi_i2,rho,maxiter,convthresh,x,lossname,dorec);
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
        %L2 = trw_bprop(model,psi_ij2,psi_i,rho,maxiter,convthresh,@loss,dorec);
        L2 = implicit_loss(model,psi_ij2,psi_i,rho,maxiter,convthresh,x,lossname,dorec);
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

drawnow
pause

fprintf('\n\n\n\n');

end

end