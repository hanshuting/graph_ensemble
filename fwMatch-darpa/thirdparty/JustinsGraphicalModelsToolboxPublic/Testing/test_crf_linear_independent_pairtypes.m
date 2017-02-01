function test_crf_linear_independent_pairtypes

ly    = 10;
lx    = 10;
nvals = 2;
nfeat = 3;
rho   = .5;
%lossname = 'ul';
loss_spec = 'trunc_ul_5';

model = gridmodel(ly,lx,nvals);

y = ceil(nvals*rand(model.nnodes,1));
x = randn(model.nnodes,nfeat);

F = randn(nvals,nfeat);
G = randn(nvals*nvals,2);

[L b_ij b_i dF dG] = crf_linear_independent_pairtypes(model,F,G,x,y,rho,loss_spec);

e = 1e-10;
for i=1:size(F,1)
    for j=1:size(F,2)
        F2 = F;
        F2(i,j) = F(i,j)+e;
        L2 = crf_linear_independent_pairtypes(model,F2,G,x,y,rho,loss_spec);
        dF2(i,j) = (1/e)*(L2-L);
    end
end

dF2
dF

for i=1:size(G,1)
    for j=1:size(G,2)
        G2 = G;
        G2(i,j) = G(i,j)+e;
        L2 = crf_linear_independent_pairtypes(model,F,G2,x,y,rho,loss_spec);
        dG2(i,j) = (1/e)*(L2-L);
    end
end

dG2
dG