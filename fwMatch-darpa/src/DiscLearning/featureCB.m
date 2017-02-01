function f = featureCB(param, x, y) 
feats = y'*param.obj.getFeatures(x);
f = sparse(feats)';
%f = param.obj.G(param.obj.beginRow(x):param.obj.endRow(x),:);
