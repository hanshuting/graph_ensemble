function yhat = constraintCB(param, model, x, trueLabel)
    theta = model.w;
   % inputFeatures = th.G(inds,:);
   %this is the margin scaling version
    
    weights = param.obj.getFeatures(x)*theta;
    weights = weights + (1- trueLabel);
    
    yhat = vecUniMatch(-weights)'; 
        
end