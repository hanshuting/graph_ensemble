% function X2 = fourier_expand(X,maxk)
% 
% C = tf_table(size(X,1),maxk+1)-1;
% 
% X2 = zeros(size(C,2),size(X,2));
% for j=1:size(C,2)
%     c = C(:,j);
%     X2(j,:) = cos(pi*c'*X);
% end

function X2 = fourier_expand(X,maxk)

if length(maxk)==2
    ksumlim = maxk(2);
    maxk = maxk(1);
end

C = tf_table(size(X,1),maxk+1)-1;

N = sum(C,1);
%max(N)
if exist('ksumlim','var')
    C = C(:,N<=ksumlim);
end

% N = sum(C~=0,1);
% C = C(:,N<=1);

X2 = zeros(2*size(C,2),size(X,2));
where = 1;
for j=1:size(C,2)
    c = C(:,j);
    X2(where,:) = cos(pi*c'*X);
    where = where+1;
    X2(where,:) = sin(pi*c'*X);
    where = where+1;
end