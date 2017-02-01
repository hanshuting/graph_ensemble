function hl = hl_factor(x);
% HL_FACTOR	Given row sums in x, computes the Huber/Law factor
%    The Minc factor is (a!)^(1/a).  The Huber/Law factor is slightly
%    larger, but asymptotically the ratio approaches 1.  Furthermore,
%    the HL factor can be used in building an acceptance/rejection
%    scheme for generating perfect matchings, while the Minc factor
%    cannot.
%
% Input Variables:
%   x            The row sum
%
% Output Variables:
%   hl           The hl factor

x = x + .5*(x==0);  % so that not taking log of 0
hl = (x > 1).*(x + .5*log(x)+exp(1)-1) + (x <= 1).*(1 + (exp(1) - 1)*x);


