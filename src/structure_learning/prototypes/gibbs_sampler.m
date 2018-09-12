% GIBBS_SAMPLER receives parameters of a Markov Random Field in the
% Ising Model and through stochastic methods samples the variables and 
% determines the marginals of each random variable.
%
% Usage
%   sampler = gibbs_sampler(node_potentials, edge_potentials, (optional) burn_in, (optional) gap)
%
% Input
%   node_potentials: column-array with the node individual potentials.
%                    Implicitly defines # of nodes = # of rows
%   edge_potentials: (# of edges)-by-3 matrix. First and second columns
%                    hold the index of the nodes conected (any order).
%                    Third column holds the edge potential.
%                    Alternatively you can pass sparse matrix.
%   burn_in:         # of samples discarded before start recording
%   gap:             # of samples discarded between recorded samples
%
% Example
%   sampler.sample(samples_to_record)
%   samples = sampler.get_samples()
%   marginals = sampler.get_marginals()

classdef gibbs_sampler < handle    
    properties (SetAccess = private, GetAccess = private)
        % CONFIGURATIONS
        node_potentials; % column array with the individual potentials of each node
        
        edge_potentials; % symmetric sparse matrix with the potentials of the edges 
                         % (activated when both nodes are 1). To avoid
                         % duplication we only set the left-bottom triangle
                         % of the matrix (e.g. column index <= line index)
        
        burn_in;    % # of samples discarded before start recording
        gap;        % # of samples discarded between recorded samples
        
        % INTERNAL VARIABLES
        N;   % node count
        samples;     % matrix of type "logical" storing the samples
        sample_count;   % samples made by now
        pseudo_loglikelihood;   % cumulative pseudo loglikelihood to check if model 
                                % is stationary before making samples (not auto
                                % but you can ask for it and set the burn_in better
                                % next time)
    end
    
    methods
        % Constructor (see docs in the beginning of class)
        function self = gibbs_sampler(node_potentials, edge_potentials, varargin);
            p = inputParser;

            % ensures it is a matrix with all nodes are between [0,N] and 
            % there are no edges from any node to itself
            valid_edges = @(x) ismatrix(x) && (max(max(x(:,1:2))) <= length(node_potentials)) && (min(min(x(:,1:2))) > 0) && all(x(:,1)-x(:,2)~=0);
            
            % convert sparse matrix to matrix of indexes and values
            if issparse(edge_potentials)
               [a b c] = find(edge_potentials);
               edge_potentials = [a b c];
            end
            
            p.addRequired('node_potentials', @isvector);
            p.addRequired('edge_potentials', valid_edges);
            p.addOptional('burn_in', 600, @isnumeric);
            p.addOptional('gap', 200, @isnumeric);
            p.parse(node_potentials, edge_potentials, varargin{:});
            
            self.node_potentials = p.Results.node_potentials;
            self.N = length(self.node_potentials);
            self.burn_in = p.Results.burn_in;
            self.gap = p.Results.gap;
            self.sample_count = 0;
            self.samples = 0;
            
            % makes sure all edges will be represented with row >= column
            aux = sort(p.Results.edge_potentials(:,1:2),2,'descend');
            edge_i = aux(:,1); 
            edge_j = aux(:,2);
            
            % create sparse matrix for edge potentials
            self.edge_potentials = sparse(edge_i, edge_j, p.Results.edge_potentials(:,3), self.N, self.N);
        end
                
        % returns full sample matrix
        function samples = get_samples(self)
            samples = self.samples;
        end
        
        % get marginals for each random variable
        function marginals = get_marginals(self)
            marginals = sum(self.samples)/self.sample_count; 
        end
        
        function pseudo_loglikelihood = get_pseudo_loglikelihood(self)
            pseudo_loglikelihood = self.pseudo_loglikelihood;
        end
        
        % start gibbs sampling
        function sample(self, samples_to_record)
            % preallocate samples matrix
            self.samples = false(samples_to_record, self.N);
            self.sample_count = 0;
            
            % iterations count involves burn in and gap
            iteration_count = self.burn_in + (samples_to_record * (1 + self.gap) - self.gap);
            
            % allocate a vector to store pseudolikelihood (we want to make
            % sure that when we record the samples they are representative
            % of the true model and the burn in has passed)
            self.pseudo_loglikelihood = zeros(1, iteration_count);
            
            % generate random state for random variables
            x = rand(1, self.N) > 0.5;
            
            for i = 1:iteration_count        
                % for each variable, fix the others and sample from the
                % distribution
                for j = 1:self.N
                    % Since all other variables are fixed, only need to worry about the node 
                    % potential of the sampled variable and edge potentials involving this variable.
                    % The ratio of p(x=1|...)/p(x=0|...) is actually a subtraction
                    % of log probabilities, but all that would be subtracted is simplified
                   
                    % init with the node potential and add all the edge potentials (where there is edge and
                    % the other variable is set). We have to do a trick because we store all values in the 
                    % lower-left triangle of the matrix
                    log_ratio = self.node_potentials(j);
                    for k=1:self.N
                        if k ~= j && x(j)
                            if k < j
                                log_ratio = log_ratio + self.edge_potentials(j,k);
                            else
                                log_ratio = log_ratio + self.edge_potentials(k,j);
                            end
                        end
                    end
                    
                    % vectorization trial that did not work:
                    %log_ratio = self.node_potentials(j) + [self.edge_potentials(j,1:j) transpose(self.edge_potentials(j+1:end,j))] * transpose(x);
                   
                    % based on the log_ratio p(x=1|...)/p(x=0|...) compute p(x=1|...)
                    prob_x_is_1 = exp(log_ratio) / (exp(log_ratio) + 1);
                   
                    % update x based on a random sample from bernoulli
                    x(j) = rand < prob_x_is_1;
                    
                    % increment pseudolikelihood for this sample
                    if x(j) == 1
                        self.pseudo_loglikelihood(i) = self.pseudo_loglikelihood(i) + log(prob_x_is_1);
                    else
                        self.pseudo_loglikelihood(i) = self.pseudo_loglikelihood(i) + log(1-prob_x_is_1);
                    end
                end
               
                % save final value of x as a sample
                if (i > self.burn_in) && (mod(i - self.burn_in - 1, self.gap + 1) == 0)
                    self.samples(self.sample_count + 1, :) = x;
                    self.sample_count = self.sample_count + 1;
                end
            end
        end
    end
    
end


