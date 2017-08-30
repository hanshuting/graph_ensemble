function [core,core_spont] = find_spont_core(mc_an,core_true)

num_stim = length(core_true);
mc_an_sz = cellfun('length',mc_an);

core = cell(num_stim,1);
min_ov = cell(num_stim,1);
for ii = 1:num_stim
    mc_ov = cellfun(@(x) length(intersect(x,core_true{ii})),mc_an);
    keep_indx = find(mc_ov>=mc_an_sz-1 & mc_an_sz>=4);
    core{ii} = unique(cell2mat(mc_an(keep_indx)));
    ov_vec = zeros(length(keep_indx));
    for jj = 1:size(ov_vec,1)
        ov_vec(jj,:) = cellfun(@(x) length(intersect(x,mc_an{keep_indx(jj)})),...
            mc_an(keep_indx));
    end
    min_ov{ii} = max(ov_vec.*double(~eye(size(ov_vec,1))),[],2);
end
ov_thresh = mean(reshape(cell2mat(min_ov),[],1));

% separable non-visual ensembles
mc_nov = cellfun(@(x) length(intersect(x,unique(cell2mat(core_true)))),mc_an);
mc_an_nov = mc_an(mc_nov<=1 & mc_an_sz>=4);
core_spont = cell(length(mc_an_nov),1);
for ii = 1:length(mc_an_nov)
    pov = cellfun(@(x) length(intersect(x,mc_an_nov{ii})),mc_an_nov);
    core_spont{ii} = unique(cell2mat(mc_an_nov(pov>=ov_thresh)));
end
core_spont = core_spont(cellfun(@(x) ~isempty(x),core_spont));

    
end