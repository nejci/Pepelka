% PUTBND: Changes the [min,max] settings for both axes to allow a buffer 
%         beyond the range of the data.  If a single matrix of point coordinates 
%         is given, columns beyond the second are ignored.
%
%     Syntax: v = putbnd(x,y,{buffer},{nocall}) 
%                   OR 
%             v = putbnd([x y],buffer,{nocall})
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

% See also putbnds() if change arguments.

% RE Strauss, 9/20/97
%   12/17/99 - allow for mixed input row/col vectors.
%   3/21/00 -  allow specification of spatial buffer size.
%   3/23/00 -  avoid initial call to axis().
%   5/23/01 -  add 'nocall' option.

function v = putbnd(x,y,buffer,nocall)
  if (nargin < 2) y = []; end;
  if (nargin < 3) buffer = []; end;
  if (nargin < 4) nocall = []; end;

  if (isempty(y) | isscalar(y))
    nocall = buffer;
    buffer = y;
    if (size(x,2)>=2)
      y = x(:,2);
      x = x(:,1);
    else
      error('  PUTBND: invalid point coordinates');
    end;
  end;

  x = x(:);
  y = y(:);
  if (length(x) ~= length(y))
    error('  PUTBND: lengths of coordinate vectors are incompatible.');
  end;

  if (isempty(buffer))
    buffer = 0.05;
  end;
  if (isempty(nocall))
    nocall = 0;
  end;

  % Remove NaN's
  indx = (isfinite(x) & isfinite(y));
  x = x(indx);
  y = y(indx);

  v = zeros(1,4);
  v(1) = min(x)-buffer*range(x);
  v(2) = max(x)+buffer*range(x);
  v(3) = min(y)-buffer*range(y);
  v(4) = max(y)+buffer*range(y);

  if (v(2)-v(1) < eps)      % No variation in x
    v(1) = x(1)-1;
    v(2) = x(1)+1;
  end;
  if (v(4)-v(3) < eps)      % No variation in y
    v(3) = y(1)-1;
    v(4) = y(1)+1;
  end;

  if (~nocall)
    axis(v);
  end;

  return;
