function T = textonify(im,C,method)

% make an image into textons

if ndims(im)==3
    im    = rgb2gray(im);
end
[ly lx] = size(im);

if strcmp(method,'MR8')
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
elseif strcmp(method,'canny');
    nclust = size(C,2);
    D      = zeros(ly,lx,nclust);
    
%     E   = double(edge(im,'canny'));
%     im2 = imresize(im,.5,'bilinear');
%     E2  = double(edge(im2,'canny'));
%     E2  = imresize(E2,size(E),'bilinear');
    
    im3 = imresize(im,.25,'bilinear');
    E3  = double(edge(im3,'canny'));
    E3  = imresize(E3,[ly lx],'nearest');
    
    %E  = reflectim(E ,50);
    %E2 = reflectim(E2,50);
    E3 = reflectim(E3,50);
    
    ED = bwdist(E3);
    
    for c=1:nclust
        who{c} = find(C(:,c));
        [who_y{c} who_x{c}] = ind2sub([11 11],who{c});
        who_y{c} = (who_y{c}-6)*4;
        who_x{c} = (who_x{c}-6)*4;
        
        delta{c} = who_y{c} + size(ED,1)*(who_x{c});
    end
    
    for c=1:nclust
        [c nclust]
        for y0=1:ly
            for x0=1:lx
%                 for i=1:length(who{c})
%                     D(y0,x0,c) = D(y0,x0,c) + ED(50+y0+who_y{c}(i),50+x0+who_x{c}(i));
%                 end
                D(y0,x0,c) = sum(ED(50+y0 + (50+x0-1)*size(ED,1) + delta{c}));
            end
        end
    end
    for c=1:nclust
        T = length(who{c})+1e-5;
        D(:,:,c) = (1+D(:,:,c)) / T;
    end
    
    [~,T] = min(D,[],3);
    keyboard
else
    error('unsupported feature type');
end