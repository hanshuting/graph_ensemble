function test_crf_linear_linear

ly    = 10;
lx    = 10;
nvals = 2;
nfeat = 3;
nedgefeat = 7;
rho   = .5;
%loss_spec = 'trunc_ul_trw_5';
%loss_spec = 'pert_ul_trw_1e-5';
loss_spec = 'pseudo';
loss_spec = 'piecewise';

model = gridmodel(ly,lx,nvals);

y = ceil(nvals*rand(model.nnodes,1));
x = randn(model.nnodes,nfeat);
z = randn(model.ncliques,nedgefeat);

F = randn(nvals  ,nfeat);
G = randn(nvals^2,nedgefeat);

[L b_ij b_i dF dG] = crf_linear_linear(model,F,G,x,z,y,rho,loss_spec);

e = 1e-10;
for i=1:size(F,1)
    for j=1:size(F,2)
        F2 = F;
        F2(i,j) = F(i,j)+e;
        L2 = crf_linear_linear(model,F2,G,x,z,y,rho,loss_spec);
        dF2(i,j) = (1/e)*(L2-L);
    end
end

dF2
dF

for i=1:size(G,1)
    for j=1:size(G,2)
        G2 = G;
        G2(i,j) = G(i,j)+e;
        L2 = crf_linear_linear(model,F,G2,x,z,y,rho,loss_spec);
        dG2(i,j) = (1/e)*(L2-L);
    end
end

dG2
dG