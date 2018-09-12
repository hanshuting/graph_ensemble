function [p] = ryser_small(A)
% RYSER_SMALL  Uses Ryser to find the permanent of matrices with small entries
%
% Input Variables
%   A		The matrix to calculate the permanent of
%
% Ouput Variables
%   p		The permanent of the matrix
%
% Author:  Jenny Law, Mark Huber
% Last Modified:  08 Aug 2007

n=length(A);
p=0;
for r=0:n-1
  % First calculate the sum over subsets of the first r rows of the matrix
  if (r==0)
    rowsum=sum(A,2);
    s=1;
    for i=1:n
      s = s * rowsum(i);
    end
  else
    % Find the next subset of columns to use
    [aa,more]=ksub_next(n,r,1:r,0);
    % Find the sum of rows excluding columns in aa
    rowsum = (sum(A,2) - sum(A(:,aa),2));
    s = prod(rowsum);
    k=1;
    while(more)
      [bb,more]=ksub_next(n,r,aa,more);
      jj1=1;
      jj2=1;
      mydiff1=0;
      mydiff2=0;
      for kk=1:r
        if(min(bb(kk)~=aa)==1)
          mydiff1(jj1)=kk;
          jj1=jj1+1;
        end
        if(min(aa(kk)~=bb)==1)
          mydiff2(jj2)=kk;
          jj2=jj2+1;
        end
      end

      mydiffsize=size(mydiff1);
      mydiffsize=max(mydiffsize);
 
      for jj=1:n
        for kk=1:mydiffsize
          rowsum(jj)=rowsum(jj)+A(jj,aa(mydiff2(kk)))-A(jj,bb(mydiff1(kk)));
        end
      end
      myprod=1;
      for jj=1:n
        myprod=myprod*rowsum(jj);
      end
      s=s+myprod;
      k=k+1;
      aa=bb;
    end
  end
  p= p+(-1)^(r)*s;
end

return