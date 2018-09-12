function perms = sample_perms(A, K)
% perms   Sample n unbiased perfect matchings (permutations) with weights A
%
%   Uses the Huber-Law algorithm (hacked their code)
%
% Input Variables:
%   A            LxL input matrix (must be nonnegative)
%   K            The number of samples (not iterations)
%
% Output Variables
%   perms        KxL matrix of K permutations.
%
% Author:  Mark Huber, www.math.duke.edu/~mhuber
% Last Modified:  February 17, 2014 by Kui Tang, http://kui-tang.com

% INITIALIZATION
n = length(A);
L = size(A, 1);

% Scale the matrix using sinkhorn balancing to 4 sig figs
[B,x,y] = sinkhorn(A,.00001);

% Create best matrix with entries in [0,1] by dividing each
% row by the largest entry
row_scale = max(B')'.^(-1);
C = diag(row_scale)*B;
save_C = C;

perms = zeros(K, L);

% Initialize counters
number_of_successes = 0;

% Generate a number of random permuations equal to iterations
while number_of_successes < K
  column = 1;
  C = save_C;
  row_sums = sum(C,2);
  
  trial_perm = zeros(1, L);
  
  while column <= n
    % find hl upper bound on permanent of C
    h = hl_factor(row_sums);
    hl = prod(h/exp(1));
    % do same thing, but remove the column under consideration;
    h2 = hl_factor(row_sums - C(:,column));
    hl2 = prod(h2/exp(1));
    % Find the probability of each row being selected
    row_probs = exp(1)*C(:,column).*hl2/hl./hl_factor(row_sums - C(:,column));
    % Select a random row to fill out the permutation
    row_pick = sum(cumsum(row_probs) < rand) + 1;
    
    if (row_pick == (n + 1))  % then failed to get a permuation
      column = n + 2;
    else
      % We did get a permutation
      % KT: save the permutation
      trial_perm(column) = row_pick;    

      % Remove current column from future consideration
      row_sums = row_sums - C(:,column);  
      % Remove row from future consideration
      C(row_pick,:) = zeros(1,n);
      column = column + 1;
    end % if
    row_sums(row_pick) = 0;
  end;
  if (column == (n + 1)) % then it was a success!        
    number_of_successes = number_of_successes + 1;    
%     assert(all(trial_perm > 0 | trial_perm <= L) , 'Success but permutation invalid!');
    perms(number_of_successes,:) = trial_perm;    
  end;
end  % for loop1
