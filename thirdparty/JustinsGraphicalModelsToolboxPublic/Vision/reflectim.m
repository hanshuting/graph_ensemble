function im2 = reflectim(im,siz)
% reflect borders of size siz
% brute force doing all 8 cases separately, since that is really faster
% than a clean solution using coordinates...

[ly lx lz] = size(im);
im2 = zeros(ly+2*siz,lx+2*siz,lz);
% for z=1:lz
% %im2(:,:,z) = zeros(ly+2*siz,lx+2*siz);
% 
% % original im
% im2(siz+1:end-siz,siz+1:end-siz,z) = im(:,:,z);
% 
% % top
% im2(1:siz,siz+1:end-siz,z) = flipud(im(1:siz,:,z));
% 
% % bot
% im2(end-siz+1:end,siz+1:end-siz,z) = flipud(im(end-siz+1:end,:,z));
% 
% % left
% im2(siz+1:end-siz,1:siz,z) = fliplr(im(:,1:siz,z));
% 
% % right
% im2(siz+1:end-siz,end-siz+1:end,z) = fliplr(im(:,end-siz+1:end,z));
% 
% % top-left
% im2(1:siz,1:siz,z) = rot90(im(1:siz,1:siz,z),2);
% 
% % top-right
% im2(1:siz,end-siz+1:end,z) = rot90(im(1:siz,end-siz+1:end,z),2);
% 
% % bot-left
% im2(end-siz+1:end,1:siz,z) = rot90(im(end-siz+1:end,1:siz,z),2);
% 
% % bot-right
% im2(end-siz+1:end,end-siz+1:end,z) = rot90(im(end-siz+1:end,end-siz+1:end,z),2);
% end
if strcmp(class(im),'double')
    for z=1:lz
        im2(:,:,z) = reflectim_helper(im(:,:,z),siz);
    end
elseif strcmp(class(im),'uint8')
    im = double(im);
    for z=1:lz
        im2(:,:,z) = reflectim_helper(im(:,:,z),siz);
    end
    im2 = uint8(im2);
else
    error('reflectim only handles doubles and uint8s (fixing this is trivial)')
end


% [X Y] = meshgrid(1-siz:lx+siz,1-siz:ly+siz);
% bads = 1;
% while bads
%     bads = find(X<1);
%     X(bads) = 1-X(bads);
% end
% bads = 1;
% while bads
%     bads = find(X>lx);
%     %X(bads) = -(X-ly)+ly
%     X(bads) = 2*lx-X(bads)+1;
% end
% bads = 1;
% while bads
%     bads = find(Y<1);
%     Y(bads) = 1-Y(bads);
% end
% bads = 1;
% while bads
%     bads = find(Y>ly);
%     Y(bads) = 2*ly-Y(bads)+1;
% end
% 
% im3 = im(Y+(X-1)*ly);
% 
% im3(1:5,1:5,:)
% im2(1:5,1:5,:)
% 
% norm(im3(:)-im2(:))