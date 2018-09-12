function plots(varargin)

if ischar(varargin{1})
    mode = varargin{1};
    varargin = {varargin{2:end}};
else
    mode = 'plot';
end

% modes
% 1 - regular
% 2 - semilogx
% 3 - semilogy
% 4 - loglog

stuff = {};
where = 1;
for i=1:length(varargin)
    p = varargin{i};
    stuff{where} = 1:length(p);
    where = where + 1;
    stuff{where} = p;
    where = where + 1;
end

if strcmp(mode,'plot')
    plot(stuff{:})
elseif strcmp(mode,'semilogx')
    semilogx(stuff{:})
elseif strcmp(mode,'semilogy')
    semilogy(stuff{:})
elseif strcmp(mode,'loglog')
    loglog(stuff{:})
end
    