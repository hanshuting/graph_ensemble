function vizBeliefs(b_i, model, saveName)
% vizBeliefs  Visualize horses from marginals

    if size(b_i, 1) == 2
        % 2 x N
        bpos = vec(b_i(2,:));
    elseif size(b_i, 2) == 2
        bpos = b_i(:,2);
    end        

    I = reshape(bpos, model.lx, model.ly);
    
    saveFile = sprintf('fig/horse_margs/%s', saveName);
    % TODO: Perhaps transpose
    figure;
    imshow(I);
    
    print('-dpdf', saveFile);
    
end