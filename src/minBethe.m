function [muMin, FBethe, FBetheHist] = minBethe(mu, A)            
% minBethe  Compute minimizer of Bethe partition function
%   [muMin, FBethe, FBetheHist] = minBethe(mu, A) where
%   
%   mu - Matrix of beliefs (for warm start)
%   A  - *Unipartite* cost matrix
%   
%   muMin  - Argmin of Bethe free energy as matrix of beliefs
%   FBethe - Minimum of Bethe free energy
%   FBetheHist - History of energy of iterates
%
    [muMin, objectives] = fwBipartite_mex(mu, A, 1e-4);
    FBetheHist = objectives;
    FBethe = FBetheHist(end);

end
