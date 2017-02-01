function F = histify(im,rez,buf)
[ly lx lz] = size(im);
im = max(im,1e-100); % make sure positive
% first of all, bin everything
bim = ceil(im*rez);
bin = ones(ly,lx);
maxbin = rez^lz;
for z=1:lz
    bin = bin + (bim(:,:,z)-1)*rez^(z-1);
end

bin2    = reflectim(bin   ,buf);

% F = zeros(ly,lx,maxbin);
% for dy=-buf:buf
%     for dx=-buf:buf
%         for y0=1:ly
%             for x0=1:lx
%                 y = buf+y0+dy;
%                 x = buf+x0+dx;
%                 b = bin2(y,x);
%                 F(y0,x0,b) = F(y0,x0,b)+1;
%             end
%         end
%     end
% end

% about twice as fast...

F  = zeros(ly,lx,maxbin);
N  = reshape(1:ly*lx,ly,lx);
y0 = 1:ly;
x0 = 1:lx;
for dy=-buf:buf
    for dx=-buf:buf
        y  = buf+y0+dy;
        x  = buf+x0+dx;
        b  = bin2(y,x);
        F(N+(ly*lx)*(b-1)) = F(N+(ly*lx)*(b-1))+1;
    end
end
