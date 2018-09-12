%% Toy Example

% Start with a model I know the answer. Markov chain:
%   x1 - x2 - x3
%   theta1 = 1
%   theta2 = 2
%   theta3 = 3
%   theta12 = -3
%   theta23 = -1
%
%   Calculating true marginals:
%   x1  x2  x3  p*Z     p
%   0   0   0   1       .0067
%   1   0   0   exp(1)  .0183
%   0   1   0   exp(2)  .0497
%   1   1   0   1       .0067
%   0   0   1   exp(3)  .1350
%   1   0   1   exp(4)  .3670
%   0   1   1   exp(4)  .3670
%   1   1   1   exp(2)  .0497
%
%   p(x1=1) = .4416
%   p(x2=1) = .4730
%   p(x3=1) = .9186

node_potentials = [1;2;3];
edge_potentials = [0 -3 0;
                   -3 0 -1;
                   0 -1 0];

sample_count = 20000;
[X, log_like] = gibbs_fast_sampler(node_potentials, edge_potentials, 3, sample_count, 10000, 200, 1);
marginals = sum(X)/sample_count
% plot(log_like)

%% Toy Example 2 (no node potentials)

% Start with a model I know the answer. Markov chain:
%   x1 - x2 - x3
%   theta1 = 0
%   theta2 = 0
%   theta3 = 0
%   theta12 = -1
%   theta23 = -1
%
%   Calculating true marginals:
%   x1  x2  x3  p*Z     
%   0   0   0   1       .170
%   1   0   0   1       .170
%   0   1   0   1       .170
%   1   1   0   exp(-1) .0626
%   0   0   1   1       .170
%   1   0   1   1       .170
%   0   1   1   exp(-1) .0626
%   1   1   1   exp(-2) .023
%
%   p(x1=1) = .4264
%   p(x2=1) = .3187
%   p(x3=1) = .4264

node_potentials = [0;0;0];
edge_potentials = [0 -1 0;
                   -1 0 -1;
                   0 -1 0];
     
sample_count = 20000;
[X, log_like] = gibbs_fast_sampler(node_potentials, edge_potentials, 3, sample_count, 10000, 200, 1);
marginals = sum(X)/sample_count
% plot(log_like)

%% Randomly Generated Graph

rng('default'); rng(2);
% generates a graph with 100 nodes and 5% edge density
node_count = 100;
[node_potentials, edge_potentials] = rand_ising(node_count, .05);

sample_count = 10000;
burn_in = 10000;
sample_interval = 100;
rand_seed = 1;
[X, log_like] = gibbs_fast_sampler(node_potentials, full(edge_potentials), node_count, sample_count, burn_in, sample_interval, rand_seed);
marginals = sum(X)/sample_count
% plot(log_like)






