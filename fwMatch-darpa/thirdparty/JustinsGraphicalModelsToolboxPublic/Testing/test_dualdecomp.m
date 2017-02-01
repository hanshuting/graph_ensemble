function test_dualdecomp

randn('state',5);

model = gridmodel(50,50,2);
rho   = .5;

% now, make some random potentials
psi_ij = 9*randn(model.nvals^2,model.ncliques);
psi_i  = randn(model.nvals  ,model.nnodes  );

fprintf('Testing that trw and dualdecomp give the same results...\n');
tic
[b_ij b_i A nreps]   = dualdecomp(model,psi_ij,psi_i,rho);
time_trw = toc;
tic
maxiter = 10000;
convthresh = 1e-6;
[b_ij2 b_i2 A2 Z_ij2 nreps2] = trw_fast(model,psi_ij,psi_i,rho,maxiter,0,convthresh);
time_trw_fast = toc;
fprintf('[Times]:\ndualdecomp: %f\ntrw_fast:   %f  \n', time_trw, time_trw_fast);
fprintf('[Difference of marginals] (should be zero) \nunivariate: %f \nbivariate:  %f \n',...
    norm(b_i(:)-b_i2(:),inf), norm(b_ij(:)-b_ij2(:),inf))
fprintf('[Approx partition functions] \ndualdecomp: %f \ntrw_fast:   %f \ndiff:       %f (should be zero) \n',...
    A,A2,A-A2);
fprintf('\n')

nreps
nreps2

end