classdef UnipartiteMatchingTestSet < handle
    properties
        M, Ys, features, opts
    end
    
    methods
        function th = UnipartiteMatchingTestSet(Ys, features, varargin)            
            th.Ys = Ys;
            th.features = features;
            th.M = length(Ys);

            p = inputParser;
            p.parse(varargin{:});
            th.opts = p.Results;
        end
        
        function err = MAPErr(th, params)    
            % wappa wappa
            err = evalThetaUnipartite(params, th.features, th.Ys);
        end
        
        function err = MPMErr(th, params)
            err = NaN;
        end
        
    end
    
end

