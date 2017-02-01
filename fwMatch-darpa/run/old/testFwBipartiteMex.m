function [dObjEnd, dObjMin, maxDBelief] = testFwBipartiteMex(T, D, varargin)    
% NOTE: Due to numerical differences, 1e-4 is probably the smallest
% tolerance under which both will behave similarly.

p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('T', @isnumeric);
p.addRequired('D', @isnumeric);

p.addParamValue('sigma', 1);
p.addParamValue('TolFun', 1e-6);
p.addParamValue('maxIter', 1000000);
p.addParamValue('compare', true);
p.addParamValue('plot', false);
p.addParamValue('isDebug', false);
p.addParamValue('plotHeadCutoff', 0.4);
p.addParamValue('saveFile', []);

p.parse(T, D, varargin{:});            
o = p.Results;    

bel    = 1/D * ones(D,D,T);
belMex = 1/D * ones(D,D,T);

A = cell(T, 1);
objs = cell(T, 1);
objsMex = cell(T, 1);

dObjEnd    = zeros(T, 1);
dObjMin    = zeros(T, 1);
maxDBelief = zeros(T, 1);

for t = 1:T    
    A{t} = o.sigma * randn(D);
    reweight = ones(D, 1);
    
    if o.compare
        [bel(:,:,t), objs{t}] = fw_permanent_simple(bel(:,:,t), A{t}, o.TolFun, o.isDebug, o.maxIter);
    end
        
    if o.isDebug
        [belMex(:,:,t), objsMex{t}] = fwBipartite_mex_debug(belMex(:,:,t), A{t}, [], o.TolFun, o.isDebug, o.maxIter);    
    else
        [belMex(:,:,t), objsMex{t}] = fwBipartite_mex_debug(belMex(:,:,t), A{t}, [], o.TolFun, o.isDebug, o.maxIter);        
    end
            
    if o.compare        
        % Cut off the first o.plotHeadCutoff of the objs series (to better
        % evaluate the tail of the series)
        startIdx = max(1, floor(o.plotHeadCutoff * min(length(objs{t}), length(objsMex{t}))));
        oIdxs = startIdx:length(objs{t});
        mIdxs = startIdx:length(objsMex{t});                

        dObjEnd(t)    = objs{t}(end) - objsMex{t}(end);
        dObjMin(t)    = min(objs{t}) - min(objsMex{t});
        maxDBelief(t) = max(abs(vec(bel(:,:,t)) - vec(belMex(:,:,t))));
        
        if o.plot            
            figure; plot(oIdxs, objs{t}(oIdxs), mIdxs, objsMex{t}(mIdxs)); legend('lab', 'mex');
            title(sprintf('FBethe; dObjEnd = %g ; dObjMin = %g', ... 
                dObjEnd(t), dObjEnd(t)));

            figure; imagesc(bel(:,:,t) - belMex(:,:,t)); colorbar;
            title(sprintf('Belief diffs; max |diff| = %g', ...
                maxDBelief(t)));
        end
    end
end

if ~isempty(o.saveFile)
    save(o.saveFile);
end

end

