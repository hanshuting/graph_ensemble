function [pred,th] = pred_from_LL(LLs,qnoise)
% make binary predictions from LL estimation with 0/1
% INPUT:
%     LLs: vector of LL_on-LL_off
%     qnoise: quantile considered as noise
% OUTPUT:
%     pred: binary vector of prediction

th = quantile(LLs(:),qnoise);
LLs_th = LLs;
LLs_th(LLs>=th) = NaN;
th = 3*nanstd(LLs_th(:))+nanmean(LLs_th(:));
if isnan(th)
    th = 0;
end
pred = LLs>th;

end