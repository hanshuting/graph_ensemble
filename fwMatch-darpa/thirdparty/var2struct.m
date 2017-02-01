function s = var2struct(varargin)
% See http://stackoverflow.com/questions/3470654/how-can-i-move-variables-into-and-out-of-a-structure-akin-to-load-and-save-in-ma
  names = arrayfun(@inputname,1:nargin,'UniformOutput',false);
  s = cell2struct(varargin,names,2);
end
