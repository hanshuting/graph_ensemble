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
edge_potentials = [1 2 -3;
                   2 3 -1];

% set seed for random number generator
rng('default'); rng(2);

sampler = gibbs_sampler(node_potentials, edge_potentials);
sampler.sample(300);
marginals = sampler.get_marginals()
%plot(sampler.get_pseudo_loglikelihood)

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
edge_potentials = [1 2 -1;
                   2 3 -1];

% set seed for random number generator
rng('default'); rng(2);

sampler = gibbs_sampler(node_potentials, edge_potentials);
sampler.sample(600);
marginals = sampler.get_marginals()
%plot(sampler.get_pseudo_loglikelihood)

%% Randomly Generated Graph

rng('default'); rng(2);
[node_potentials, edge_potentials] = rand_ising(100, .05);
sampler = gibbs_sampler(node_potentials, edge_potentials);
sampler.sample(300); % takes ~3 minutes in regular laptop
marginals = sampler.get_marginals()
%plot(sampler.get_pseudo_loglikelihood)

