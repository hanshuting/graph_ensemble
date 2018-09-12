function featurizeImages(nvals, featParams, edgeParams, imageFiles, labelFiles, saveFile)

N = length(imageFiles);

yNodes = cell(N, 1);
yEdges = cell(N, 1);
Ns     = zeros(N, 1);
feats  = cell(N, 1);
edges  = cell(N, 1);
efeats = cell(N, 1);

for n = 1:N
    %% One sample        
    fprintf('Featurizing %d / %d\n', n, N);

    % load data
    % Because the labels are stored in jpeg, a lossy compression format, we
    % may have some noise, so we threshold above the midpoint. This
    % unfortunately means our ground truth is noisy. Go figure...
    
    % Justin's convention.
    labels0 = double(imread(labelFiles{n}) == 2);
    im      = double(imread(imageFiles{n})) / 255;

    % In case we get any grayscale images, convert to RGB (just repeat).
    if ndims(im) == 2
        im = repmat(im, 1, 1, 3);
    end
    
    % featurize
    feats{n}  = featurize_im(im,featParams);
    % reshape features
    [ly, lx, lz] = size(feats{n});
    feats{n}  = reshape(feats{n},ly*lx,lz);
    
    % make the grid model (we can't hash, since different sizes)
    model = gridmodel(ly,lx,nvals);    
    edges{n}  = int32(model.pairs');
    
    [yNodes{n}, yEdges{n}] = overcompleteLabels(labels0(:) + 1, 2, edges{n});        
    
    efeats{n} = edgeify_im(im,edgeParams,model.pairs,model.pairtype);                        
    
    Ns(n) = model.nnodes;
end

Ut = single(cell2mat(feats))';
Vt = single(cell2mat(efeats))';
YN = single(cell2mat(yNodes));
YE = single(cell2mat(yEdges));

save(saveFile, 'Ut', 'Vt', 'YN', 'YE', 'Ns', 'edges');
end


