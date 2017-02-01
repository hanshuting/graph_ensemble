function T = MR8textonify(im,C)

% make an image into textons

if ndims(im)==3
    im    = rgb2gray(im);
end
[ly lx] = size(im);

f = MR8fast(im);

% distance (squared) of each pixel to each cluster
D = zeros(8,size(f,2));
for i=1:size(C,2)
    c = C(:,i);
    c = repmat(c,1,size(f,2));
    D(i,:) = sum((c-f).^2,1);
end

[~,T] = min(D,[],1);
T = reshape(T,[ly lx]);
