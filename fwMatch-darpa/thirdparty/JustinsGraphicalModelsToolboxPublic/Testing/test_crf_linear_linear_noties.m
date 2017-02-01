function test_crf_linear_linear_noties

ly    = 5;
lx    = 5;
nvals = 2;
nfeat = 3;
nedgefeat = 7;
rho   = .5;
loss_spec = 'trunc_ul_trw_5';

model = gridmodel(ly,lx,nvals);

y = ceil(nvals*rand(model.nnodes,1));
x = randn(model.nnodes,nfeat);
z = randn(model.ncliques,nedgefeat);

F = randn(nvals  ,nfeat    ,model.nnodes);
G = randn(nvals^2,nedgefeat,model.ncliques);

[L b_ij b_i dF dG] = crf_linear_linear_noties(model,F,G,x,z,y,rho,loss_spec);

e = 1e-7;
for i=1:size(F,1)
    for j=1:size(F,2)
        for k=1:size(F,3)
            F2 = F;
            F2(i,j,k) = F(i,j,k)+e;
            L2 = crf_linear_linear_noties(model,F2,G,x,z,y,rho,loss_spec);
            dF2(i,j,k) = (1/e)*(L2-L);
        end
    end
end

%dF2
%dF
'difference of dF and numerical:'
norm(dF(:)-dF2(:))

for i=1:size(G,1)
    for j=1:size(G,2)
        for k=1:size(G,3)
            G2 = G;
            G2(i,j,k) = G(i,j,k)+e;
            L2 = crf_linear_linear_noties(model,F,G2,x,z,y,rho,loss_spec);
            dG2(i,j,k) = (1/e)*(L2-L);
        end
    end
end

%dG2
%dG
'difference of dG and numerical:'
norm(dG(:)-dG2(:))
