function [ output_args ] = evalBipartite(theta, Y, features)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
function [learnLoss, noLearnLoss] = evalGraphMatch(theta, st)
% evalGraphMatch  Compute normalized Hamming loss (percent wrong match)
%
%   st contains Fsquare and FnegExpNegSquare
    
    [M, K] = size(st.Fsquare);    
    D = size(st.Fsquare{1,1}, 1);
    
    % In the houses data, the ground truth is the identity permutation:
    % "Each frame in this sequence has been hand-labeled, with the same 30
    %  landmarks identified in each frame"    
    Ytrue = 1:D;    

    % TODO: Switcharoo all over to CSA
    lapjvTol = 1e-4;
    
    learnErrs   = 0;    
    noLearnErrs = 0;    
    for m = 1:M
        % Construct specific A matrix for sample m        
        Alearn   = zeros(D, D);            
        AnoLearn = zeros(D, D);

        % TODO: Possible sign error
        for k = 1:K            
            Alearn   = Alearn - theta(k) .* st.Fsquare{m,k};
            AnoLearn = AnoLearn + st.FnegExpNegSquare{m,k};
        end                

        YhatLearn   = csaAssignPerm(Alearn');
        YhatNoLearn = csaAssignPerm(AnoLearn');

        fprintf('Learn   [%d]: %s\n', m, num2str(YhatLearn));
        fprintf('NoLearn [%d]: %s\n', m, num2str(YhatNoLearn));
        
        learnErrs   = learnErrs   + sum(Ytrue ~= YhatLearn);
        noLearnErrs = noLearnErrs + sum(Ytrue ~= YhatNoLearn);
    end

    learnLoss   = learnErrs   / (D * M);
    noLearnLoss = noLearnErrs / (D * M);
end



end

