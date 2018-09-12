function o = log_sum_exp(A,dim)
% usage:
%
% o = log_sum_exp(A,dim)
%
% so you have a matrix A
% and you want to compute
% o = log(sum(exp(A),dim))
% but you want to to it with out exp blowing outside the range of numbers
% ieee floating point can represent.
%
% author justin domke, email (author's last name)@cs.umd.edu

if nargin < 2
    dim = 1;
end

C = -max(A,[],dim);

sizes      = 1+0*size(A);
sizes(dim) = size(A,dim);

C2 = repmat(C,sizes);
o  = log(sum(exp(A+C2),dim)) - C;
