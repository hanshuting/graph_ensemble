function [err Lall] = eval_crf_mpm_all(p,feats,efeats,labels,models,loss_spec,crf_type,rho)

N = length(feats);
mistakes = 0;
points   = 0;
Lall     = 0;
for n = 1:N
    [b_i ~ L] = eval_crf(p,feats{n},efeats{n},models{n},loss_spec,crf_type,rho);

    % b_i is nvals x nnodes.

    [~, pred] = max(b_i, [], 1);
    pred     = pred.';
    mistakes = mistakes + sum(pred ~= labels{n});
    points   = points   + length(labels{n});
    Lall     = Lall + L;
end

err = mistakes / points;

end
