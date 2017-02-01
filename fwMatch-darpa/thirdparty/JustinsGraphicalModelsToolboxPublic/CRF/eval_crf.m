function [b_i b_ij L] = eval_crf(p,feat,efeat,model,loss_spec,crf_type,rho)

if isempty(efeat)
    npairs = size(model.pairs,1);
    efeat  = ones(npairs,1);
end

x = feat;
y = zeros(model.nnodes,1);
z = efeat;
        
if strcmp(crf_type,'linear_linear')
    [L, b_ij b_i] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,loss_spec);
else
    % trivial to add others here...
    error('unsupported crf type')
end