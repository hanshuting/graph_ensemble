function per_estimate = permanent_estimate(A,iterations);
% PERMANENT_ESTIMATE   Estimates the permanent by Monte Carlo methods
%    Generates an estimate of the permanent of a matrix by 
%    using the Huber/Law algorithm
%
% Input Variables:
%   A            The input matrix (must be nonnegative)
%   iterations   The number of iterations to use generating the estimate
%
% Output Variables
%   per_estimate	The estimate of the permanent
%
% Author:  Mark Huber, www.math.duke.edu/~mhuber
% Last Modified:  August 07, 2007

% INITIALIZATION
n = length(A);

% Scale the matrix using sinkhorn balancing to 4 sig figs
[B,x,y] = sinkhorn(A,.00001);

% Create best matrix with entries in [0,1] by dividing each
% row by the largest entry
row_scale = max(B')'.^(-1);
C = diag(row_scale)*B;
save_C = C;

% Initialize counters
number_of_successes = 0;

% Generate a number of random permuations equal to iterations
for loop1=1:iterations
  column = 1;
  C = save_C;
  row_sums = sum(C,2);
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
    row_pick = sum(cumsum(row_probs) < unifrnd(0,1)) + 1;
    if (row_pick == (n + 1))  % then failed to get a permuation
      column = n + 2;
    else
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
  end;
end  % for loop1

disp(sprintf('Number_of_successes: %d',number_of_successes));
% Get HL upper bound on permanent of C
C = save_C;
row_sums = sum(C,2);
hl_C = prod(hl_factor(row_sums)/exp(1));
per_estimate = hl_C*number_of_successes/iterations;
per_estimate = per_estimate/prod(row_scale)/prod(x)/prod(y);
