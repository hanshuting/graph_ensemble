classdef SparseCtor
    % SparseCtor  Imperatively construct a sparse matrix.
    
    properties
        ri, rj, rw, size, capacity, rows, cols
    end
    
    methods
        function th = SparseCtor(varargin)
            p = inputParser;
            % DO NOT USE CAPACITY FOR NOW -- just directly index the
            % out-of-range element to create it.
            p.addParamValue('capacity', -1);
            p.addParamValue('rows', []);
            p.addParamValue('cols', []);
            
            p.parse(varargin{:});
            
            
            th.capacity = p.Results.capacity;
            th.rows = p.Results.rows;
            th.cols = p.Results.cols;            
            
            th.size = 0;
            th.ri = [];
            th.rj = [];
            th.rw = [];            
        end
        
        function th = realloc(th)
        end
                
        function th = subsasgn(th, idxs, val)
            if ~strcmp('()', idxs.type)
                error('SparseCtor:subsasgn', 'Only array () indexing supported.');
            end
            
            % TODO: Implement the full syntax later, e.g. allow
            % A([1 2 3], [3 2]) = reshape(1:6, 3, 2) to do the right thing.
            
            [is, js] = idxs.subs{1:2};
            
            if length(is) > 1 || length(js) > 1
                error('SparseCtor:subsasgn', 'Only singleton indexing supported now.');
            end
            
%             if th.size == th.capacity
%                 th.realloc();
%             end
            
            % Now we are guaranteed to have free space
            th.size = th.size + 1;
            
            th.ri(th.size) = is(1);
            th.rj(th.size) = js(1);
            th.rw(th.size) = val;
        end        
        
        function S = makeSparse(th, varargin)
            p = inputParser;
            p.addParamValue('symmetrize', false);
            
            p.parse(varargin{:});
            
            if p.Results.symmetrize
            end
            
            % TODO: Implement doubling realloc.
            if ~isempty(th.rows) && ~isempty(th.cols)
                S = sparse(th.ri, th.rj, th.rw, th.rows, th.cols);
            else
                S = sparse(th.ri, th.rj, th.rw);
            end            
        end        
    end
    
end

