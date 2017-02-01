function ll = treegmLogprobBinary(model, X)
% log probabiltiy of a fully observed discrete data vector under a tree model
% LL(n) = log p(X(n,:) | params)

% This file is from pmtk3.googlecode.com

CPDs = model.CPDs;
%G = model.G;
[N d] = size(X);

ll = zeros(N,1);
% We avoid iterating over data cases
for i=1:d
   %j = parents(G, i);
   j = model.pa(i);
   CPT = CPDs{i};
   if j==0 % no parent
      ll = ll + (1-X(:,i))*log(CPT(1)+eps) + X(:,i)*log(CPT(2)+eps);
%           - logZ(CPT); % - logZ
   else
      ll = ll + (1-X(:,j)).*(1-X(:,i))*log(CPT(1,1)+eps) ...
          + (1-X(:,j)).*X(:,i)*log(CPT(1,2)+eps) ...
          + X(:,j).*(1-X(:,i))*log(CPT(2,1)+eps) ...
          + X(:,j).*X(:,i)*log(CPT(2,2)+eps);
%           - logZ(CPT); % - logZ
   end
end

end

function lz = logZ(P)
elz = sum(P(:));
if length(P(:)) == 2
    elz = elz + sqrt(P(1)*P(2));
else
    elz = elz + sqrt(P(1)*P(2)) + sqrt(P(1)*P(3)) ...
        + sqrt(P(4)*P(2)) + sqrt(P(4)*P(3));
    elz = elz + sqrt(sqrt(P(1)*P(2)*P(3)*P(4)));
end
lz = log(elz + eps);
end