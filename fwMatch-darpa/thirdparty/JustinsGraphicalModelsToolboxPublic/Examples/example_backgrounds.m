function example_backgrounds

imsdir = '~/Datasets/iccv09Data/images/'; % Change this to fit your system!
labdir = '~/Datasets/iccv09Data/labels/'; % Change this to fit your system!
nvals  = 8;
rez    = .2; % how much to reduce resolution
rho    = .5; % (1 = loopy belief propagation) (.5 = tree-reweighted belief propagation)

ims_names = dir([imsdir '*.jpg']);
lab_names = dir([labdir '*regions.txt']);

% labels are listed as an array of integers 0-7 with negative for unlabeled in a text file
% must convert to our representation of 1-8 with 0 for unlabeled

N = length(ims_names);
ims    = cell(N,1);
labels = cell(N,1);

feat_params = {{'patches',0},{'position',1},{'fourier',1},{'hog',8}};

fprintf('loading data and computing feature maps...\n');
parfor n=1:N
% for n=1:N
    fprintf('Featurizing %d / %d\n', n, N);
    % load data
    lab = importdata([labdir lab_names(n).name]);
    im  = double(imread(([imsdir ims_names(n).name])))/255;
    ims{n}  = im;
    labels0{n} = max(0,lab+1);
    
    % compute features
    feats{n}  = featurize_im(ims{n},feat_params);
    
    % reduce resolution for speed
    ims{n}    = imresize(ims{n}   ,rez,'bilinear');
    feats{n}  = imresize(feats{n} ,rez,'bilinear');
    labels{n} = imresize(labels0{n},rez,'nearest');
    
    % reshape features
    [ly lx lz] = size(feats{n});
    feats{n} = reshape(feats{n},ly*lx,lz);
end 

% the images come in slightly different sizes, so we need to make many models
% use a "hashing" strategy to not rebuild.  Start with empty giant array
model_hash = repmat({[]},1000,1000);
fprintf('building models...\n')
for n=1:N
    [ly lx lz] = size(ims{n});
    if isempty(model_hash{ly,lx});
        model_hash{ly,lx} = gridmodel(ly,lx,nvals);
    end
end
models = cell(N,1);
for n=1:N
    [ly lx lz] = size(ims{n});
    models{n} = model_hash{ly,lx};
end

fprintf('computing edge features...\n')
edge_params = {{'const'},{'diffthresh'},{'pairtypes'}};
efeats = cell(N,1);
parfor n=1:N
    efeats{n} = edgeify_im(ims{n},edge_params,models{n}.pairs,models{n}.pairtype);
end

fprintf('splitting data into a training and a test set...\n')
% split everything into a training and test set

k = 1;
[who_train who_test] = kfold_sets(N,5,k);

ims_train     = ims(who_train);
feats_train   = feats(who_train);
efeats_train  = efeats(who_train);
labels_train  = labels(who_train);
labels0_train = labels0(who_train);
models_train  = models(who_train);

ims_test     = ims(who_test);
feats_test   = feats(who_test);
efeats_test  = efeats(who_test);
labels_test  = labels(who_test);
labels0_test = labels0(who_test);
models_test  = models(who_test);


    % visualization function
    function viz(b_i)
        % here, b_i is a cell array of size nvals x nvars
        M = 5;
        for n=1:M
            [ly lx lz] = size(ims_train{n});
            subplot(3,M,n    ); miximshow(reshape(b_i{n}',ly,lx,nvals),nvals);
            subplot(3,M,n+  M); imshow(ims_train{n})
            subplot(3,M,n+2*M); miximshow(reshape(labels_train{n},ly,lx),nvals);
            
        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end
    
% SAVE PRECOMPUTED FEATURES BEFORE HERE AND COPY-PASTE THE PART BELOW.
save -v7.3 example_backgrounds_featurize


fprintf('training the model (this is slow!)...\n')
loss_spec = 'trunc_cl_trwpll_5';
crf_type  = 'linear_linear';
options.viz         = @viz;
options.print_times = 0; % since this is so slow, print stuff to screen
options.gradual     = 1; % use gradual fitting
options.maxiter     = 1000;
options.rho         = rho;
options.reg         = 1e-4;
options.opt_display = 0;
p = train_crf(feats_train,efeats_train,labels_train,models_train,loss_spec,crf_type,options);
save p p
%p = [];
%load p

% use this to train using all data 
%p = train_crf(feats,efeats,labels,models,loss_spec,crf_type,options);


fprintf('get the marginals for test images...\n');
close all
for n=1:length(feats_test)
    [b_i b_ij] = eval_crf(p,feats_test{n},efeats_test{n},models_test{n},loss_spec,crf_type,rho);
    
    [ly lx lz] = size(labels_test{n});
    [~,x_pred] = max(b_i,[],1);
    x_pred = reshape(x_pred,ly,lx);
    
    [ly lx lz] = size(labels0_test{n});
    x       = labels0_test{n};
    % upsample predicted images to full resolution
    x_pred  = imresize(x_pred,size(x),'nearest');
    E(n)   = sum(x_pred(x(:)>0)~=x(x(:)>0));
    T(n)   = sum(x(:)>0);
    
    [ly lx lz] = size(ims_test{n});
    subplot(2,3,1)
    miximshow(reshape(b_i',ly,lx,nvals),nvals);
    subplot(2,3,2)
    imshow(ims_test{n})
    subplot(2,3,3)
    miximshow(reshape(labels_test{n},ly,lx),nvals);
    
    [ly lx lz] = size(labels0_test{n});
    subplot(2,3,4)
    miximshow(reshape(x_pred,ly,lx),nvals);
    subplot(2,3,5)
    imshow(ims_test{n})
    subplot(2,3,6)
    miximshow(reshape(labels0_test{n},ly,lx),nvals);
    drawnow
end
fprintf('total pixelwise error on test data: %f \n', sum(E)/sum(T))

% colormap-- taken from http://dags.stanford.edu/projects/scenedataset.html
cmap = [.5  .5  .5;
        .5  .5   0;
        .5  .25 .5;
        0   .5   0;
        0    0   .5;
        .5   0   0;
        .5  .31  0;
        1   .5   0];
        
% takes about .62s to compute features and 
% .72s to compue marginals
    
% for fun, run on a video
clf
vid_read  = VideoReader('~/Datasets/georgetown_driving.mp4');
vid_write = VideoWriter('marginals'); open(vid_write);
for i=1:3:vid_read.NumberOfFrames/2,
    % load data, make features, etc.
    im0   = double(read(vid_read,i))/255;
    feats = featurize_im(im0,feat_params);
    im    = imresize(im0,.2,'bilinear');
    feats = imresize(feats,.2,'bilinear');
    [ly lx lz] = size(feats);
    feats = reshape(feats,ly*lx,lz);
    if i==1
        [ly lx lz] = size(im);
        model = gridmodel(ly,lx,nvals);
    end
    efeats = edgeify_im(im,edge_params,model.pairs,model.pairtype);
    
    % do inference
    [b_i b_ij] = eval_crf(p,feats,efeats,model,loss_spec,crf_type,rho);
    
    % interpolate beliefs to original image, get predictions
    [ly0 lx0 lz] = size(im0);
    b_i = reshape(b_i',[ly lx nvals]);
    b_i = imresize(b_i,[ly0 lx0],'bilinear');
    [~,x_pred] = max(b_i,[],3);

    %preds = miximshow(x_pred,nvals);
    im_gray = repmat(rgb2gray(im0),[1 1 3]);
    
    colormap(cmap)
    preds = miximshow(b_i,nvals);
    im_mix = .25*im_gray + .75*preds;
    imshow(im_mix)
    colorbar('Location','South','XTickLabel',...
    {'sky','tree','road','grass','water','bldg','mntn','fg obj'});
    
    % write output video
    currFrame = getframe;
    writeVideo(vid_write,currFrame);

    drawnow;
end
close(vid_write);

end
