function [B,x,y] = sinkhorn(A,epsilon);
% SINKHORN     Sinkhorn balances a matrix
%    In sinkhorn balancing, the row sums become one by scaling the rows.
%    Then the columns are scaled to one, then the rows, and so one.  This
%    process continues until both the rows and columns are within epsilon 
%    of 1.  This process is guaranteed to terminate after a polynomial 
%    number of iterations.
%
% Input Variables:
%   A          The matrix to be scaled
%   epsilon    Ends scaling when row sums are within epsilon of 1
%
% Output Variables
%   B		The scaled matrix
%   x          The scaling for the rows
%   y          The scaling for the columns
%
% Author:  Mark Huber, www.math.duke.edu/~mhuber
% Last Modified:  August 07, 2007

% INITIALIZATION
B = A;
n = length(A);
x = ones(n,1);
y = x';

% Start scaling until row sums and column sums are within epsilon of 1.
row_sum = sum(A,2);
while (max(abs(row_sum' - 1)) > epsilon)
  x = x .* row_sum.^(-1);
  B = diag(row_sum.^(-1))*B;
  col_sum = sum(B,1);
  y = y .* col_sum.^(-1);
  B = B*diag(col_sum.^(-1));
  row_sum = sum(B,2);
end % while loop
