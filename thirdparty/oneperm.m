function perm = oneperm(N,M)
% ONEPERM - obtain the M-th permutation of (1:N)
%   perm = ONEPERM(N,M) returns the m-th permutation of the sorted list
%   of all permutations from PERMS, where M=1 corresponds to identity 
%   permutation. N, M are non-negative scalar, perm has size 1-by-N.
%
%   See also PERMS
%        and NPERMUTEK, RECPERMS, NEXTPERM, PERMS1 on the File Exchange
%
% Algorithm: For given N and M, where 1 <= M <= N!, the M-th 
%            permutation of N objects is closely related to the 
%            factoradic of M; see factoradic on the File Exchange.
%            To convert the factoradic into a permutation follow these
%            steps
%            
%            For decreasing i
%                If element(j)>=element(i) ; where j>i
%                    element(j) increase by one.
%
%            The result will be a permutation of (0:N-1).
%            Add 1 to yield the permutation of (1:N).

% for Matlab (should work for most versions)
% version 1.0 (Feb 2009)
% (c) Darren Rowland
% email: darrenjrowland@hotmail.com
%
% Keywords: single permutation

error(nargchk(2,2,nargin));
if numel(N) ~= 1 || N <= 0 || N ~= round(N)
  error('oneperm:InvalidArg1',...
        'The first input has to be a non-negative integer');
end
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
  error('oneperm:InvalidArg2',...
        'The second input has to be a non-negative integer');
end

if M>factorial(N)
     error('oneperm:largeM','M should not exceed N!');
end
% convert M to zero-based
M = M-1;
perm = factoradic(M,N);

for ii= N-1:-1:1
%     for jj = ii+1:N
%         if(perm(jj)>=perm(ii))
%             perm(jj) = perm(jj) + 1;
%         end
%     end
    perm(ii+1:N) = perm(ii+1:N) + (perm(ii+1:N)>=perm(ii));
end
% convert permutation to one-based (from zero-based)
perm = perm + 1;

% LOCAL FUNCTIONS

function F = factoradic(M,N)
% See FACTORADIC, available on the File Exchange
F = zeros(1,N);
jj = 2;
while M~=0
    F(N-jj+1) = mod(M,jj);
    M = floor(M/jj);
    jj = jj+1;
end