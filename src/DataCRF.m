classdef DataCRF
% Converts structure + samples in the overcomplete parametrization
% necessary to run the parameter estimation code
    
    properties 
        Ns, edges, YE, YN, Ut, Vt
    end
    
    methods
        function obj = DataCRF(graphs,samples)
            
            % Initialization
            Ns = [];
            edges = {};
            YE = [];
            YN = [];
            Ut = [];
            Vt = [];
            
            for i = 1:length(graphs)
                graph = graphs{i};
                X = samples{i};
                X = (X == 1) + 1;
                
                M = size(X,1);
                N = size(X,2);
                
                obj.Ns = [obj.Ns;repmat(N,M,1)];
                
                % Get the edges from the graph, such that for (i,j) \in
                % edges, i <j
                [r,c] = find(graph);
                edgeList = [c';r'];
                [r,c] = find(edgeList(1,:) >= edgeList(2,:));
                edgeList(:,c) = [];
                obj.edges = [obj.edges;repmat({edgeList},M,1)];
                
                % E is number of edges
                E = size(edgeList,2);
                
                % YN and YE
                for m=1:M
                    [YNodes,YEdges] = overcompleteLabels(X(m,:)',2,edgeList);
                    obj.YN = [obj.YN;YNodes];
                    obj.YE = [obj.YE;YEdges];
                end
                
                % Ut and Vt 
                obj.Ut = eye(N); % D * N where D = N (D = #features) 
                obj.Vt = zeros(N*(N-1)/2,E); % D * E where E is number of edges and D = N*(N-1)/2 (all possible edges)
                for e=1:E
                    obj.Vt(obj.getEdgeId(edgeList(:,e),N),e) = 1;
                end
                obj.Ut = sparse(obj.Ut);
                obj.Vt = sparse(obj.Vt);
                
                % Stack Ut and Vt
                obj.Ut = repmat(obj.Ut,1,M);
                obj.Vt = repmat(obj.Vt,1,M);
            end
        end
        
        function id = getEdgeId(obj,edge,N)
            i = edge(1);
            j = edge(2);
            id = N*(i-1) - sum(1:(i-1)) + (j-i);
        end
    end
end


