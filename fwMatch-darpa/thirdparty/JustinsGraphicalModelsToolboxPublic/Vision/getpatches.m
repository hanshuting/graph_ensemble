function [X Y] = getpatches(imsdir,labelsdir,feat_params,npatches,maxims)

%function [params times fvals] = train_labelers(loss_spec,params0,testmode,...
%    imsdir,labelsdir,nvals,feat_params,edge_params,maxims)

Aims    = [dir([imsdir    '/*.png']);dir([imsdir    '/*.jpg']);dir([imsdir    '/*.bmp'])];
Alabels = [dir([labelsdir '/*.png']);dir([labelsdir '/*.jpg']);dir([labelsdir '/*.bmp']);];

if exist('maxims','var')
    Aims    = Aims   (1:maxims);
    Alabels = Alabels(1:maxims);
end

where = 0;
for n=1:length(Aims)
    im    = imread([imsdir    '/' Aims(n).name]);
    im    = double(im)/255;
    if size(im,3)==1, im=repmat(im,[1 1 3]); end;
    label  = double(imread([labelsdir '/' Alabels(n).name]));
    [ly lx] = size(label);
    %model   = modelAR{ly,lx};
    feats   = featurize_im(im,feat_params);
    
    if n==1
        Y = zeros(1,npatches);
        X = zeros(size(feats,3),npatches);
    end
    
    patches_per_im = ceil((npatches-where)/(length(Aims)-n+1));
    for m=1:patches_per_im
        y = 0;
        where = where+1;
        while y==0
            y_cor = ceil(ly*rand);
            x_cor = ceil(lx*rand);
            y = label(y_cor,x_cor);
            x = feats(y_cor,x_cor,:);
            Y(where)   = y;
            X(:,where) = x;
        end
    end
end