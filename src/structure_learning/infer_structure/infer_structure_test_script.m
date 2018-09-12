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
sampler.sample(3000);
samples = sampler.get_samples();
graph = infer_structure(samples, 0.1);

%% Randomly Generated Graph

% rng('default'); rng(2);
% [node_potentials edge_potentials] = randIsing(100, .05);
% sampler.sample(3000); % takes ~3 minutes in regular laptop
% marginals = sampler.get_marginals()


