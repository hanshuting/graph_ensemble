function [h, cornerA, cornerB] = plotImageHelper(cornerFileA, imageFileA, cornerFileB, imageFileB)
    
    imgA = imread(imageFileA, 'png');
    imgB = imread(imageFileB, 'png');
    
    imgAB = double([imgA imgB]);
    % Normalize to make the rest easier
    imgAB = imgAB ./ max(imgAB(:));
    
    BxOffset = size(imgA, 2);
    
    cornerA = importdata(cornerFileA);    
    cornerB = importdata(cornerFileB);
    cornerB(:,1) = cornerB(:,1) + BxOffset;
    
    h = figure;
    RGB = cat(3,imgAB,imgAB,imgAB);
    imshow(RGB);
    axis off;
    hold on;

end
