% PUTBNDS: Changes the [min,max] settings for both axes to allow a buffer 
%          beyond the range of the data.  If a single matrix of point coordinates 
%          is given, columns beyond the second are ignored.
%
%     Syntax: v = putbnds(x,y,{buffer},{nocall}) 
%                   OR 
%             v = putbnds([x y],buffer,{nocall})
%
%          x =      vector of x coordinates.
%          y =      vector of y coordinates.
%          buffer = optional buffer size, as proportion of ranges 
%                     [default = 0.05].
%          nocall = optional boolean flag indicating, if true, that the axis 
%                     settings are to be returned but that the current plot is 
%                     to be left unaltered [default = 0].
%          -------------------------------------------------------------------
%          v =      axis bounds: [xmin xmax ymin ymax].
%

% RE Strauss, 2/26/00
%   5/23/01 - change arguments to match change in putbnd().

function v = putbnds(x,y,buffer,nocall)
  if (nargin < 2) y = []; end;
  if (nargin < 3) buffer = []; end;
  if (nargin < 4) nocall = []; end;

  v = putbnd(x,y,buffer,nocall);

  return;
