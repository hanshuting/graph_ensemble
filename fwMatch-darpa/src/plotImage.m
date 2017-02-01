function h = plotImageMatch(theta, pairTxt)

    % corner feature image adj
    
    fpair = fopen(pairTxt);
    files = textscan(fpair, '%s %s %s %s', inf); 
    fclose(fpair);
% 
%     
%     t = 1;
%     while ~feof(fpair)
%         for i = 1:4
%             corner  = fscanf(fpair, '%s ');
%             feature = fscanf(fpair, '%s ');
%             image   = fscanf(fpair, '%s ');
%             adj     = fscanf(fpair, '%s ');
%         end
%         files{t} = var2struct(corner, feature, image, adj);
%         t = t + 1;
%     end
% 
%     files = cell2mat(files)
        
    files
    

end

