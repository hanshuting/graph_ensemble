function mixim=miximshow(probs,num_labels,f)

[ly lx lz] = size(probs);

if lz==1
    probs0 = probs;
    probs = zeros(ly,lx,num_labels);
    for l=1:num_labels
        probs(:,:,l)=(probs0==l);
    end
    lz = num_labels;
end

%size(probs)
%num_labels
%size(f)
%nargin

if nargin==3
    %probs = exp(f);
    %Z = sum(probs,3);
    %for z=1:lz
    %    probs(:,:,z) = probs(:,:,z)./Z;
    %end
    probs = exp(f-repmat(log_sum_exp(f,3),[1 1 lz]));    
end

[ly lx lz] = size(probs);

% if lz==2
%     mixim = repmat(probs(:,:,2),[1 1 3]);
% elseif lz==5    
%     cmap   = [[.9 0 0]; [.4 .4 .4]; [0 0 1]; [0 1 0]; [.6 .3 1]];
%     layers = {'building', 'road', 'sky', 'foliage', 'cars'};
% 
%     % display mixed
%     mixim = zeros(ly,lx,3);
%     for i=1:lz
%         for c=1:3
%             mixim(:,:,c) = mixim(:,:,c) + cmap(i,c)*probs(:,:,i);
%         end
%     end
% else
%     map = colormap('hsv');
%     l = size(map,1);
%     mixim = zeros(ly,lx,3);
%     for i=1:lz
%         for c=1:3
%             mixim(:,:,c) = mixim(:,:,c) + map(ceil(i*l/lz),c)*probs(:,:,i);
%         end
%     end
%         
% end

% compute mean
% mixim = zeros(ly,lx);
% for z=1:lz
%     mixim = mixim + z*probs(:,:,z);
% end

mixim = zeros(ly,lx,3);
for z=1:lz
    %colormap gray
    map = colormap;
    %col = map(ceil(z/lz*size(map,1)),:); % old way-- tried and true
    col = map(max(1,ceil((z-1)/(lz-1)*size(map,1))),:);
    for chan=1:3
        mixim(:,:,chan) = mixim(:,:,chan) + col(chan)*probs(:,:,z);
    end
end

mixim = max(mixim,0);
mixim = min(mixim,1);

if nargout==0
    %imagesc(mixim,[1 num_labels]);
    %axis image
    %axis off
    imshow(mixim);
end