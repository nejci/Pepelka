function [f,dfdx,dfdy]=surf_ellip(np,d,hn,ix,iy,res)
% [f,dfdx,dfdy]=surf_ellip(np,d,hn,ix,iy,res)
%
% surf_ellip computes the surface function and partial
% derivatives of a modified gaussian, where the variance
% is influenced by an elliptic arc in the corresponding quadrant.
%
% Inputs:
% np - (vector) neuron position of central neuron, (1 x 2) vector
% d - (matrix) normalized distance matrix
% hn - (double) normalized frequency of neuron np
% ix - (vector) x dimension of euclidian plane, [x_low,x_high]
% iy - (vector) y dimension of euclidian plane, [y_low,y_high]
% res - (double) resolution used for surface 
%
% Outputs:
% sf - (matrix) surface function
% dfdx, dfdy - (matrix) partial derivatives
% 
% $Id: surf_ellip.m 481 2005-03-24 08:58:38Z dome $
% D. Brugger, 15 March 2005
% algo/surf_ellip.m

% prepare grids
[X,Y]=meshgrid(ix(1):res:ix(2),iy(1):res:iy(2));
f=zeros(size(X)); 
dfdx=zeros(size(X)); dfdy=zeros(size(X));
iX=np(1)-X; jY=np(2)-Y;
iX2=iX.^2; jY2=jY.^2;
d2=(1-d).^2;
pos=mesh_idx(np,res);
dim=mesh_idx([ix(2) iy(2)],res);
x_h=dim(1); y_h=dim(2);

if(~(on_right_border(np,ix) || on_upper_border(np,iy)))
  % compute first quadrant
%  fprintf('computing first quadrant...\n')
  [yrange,xrange]=quadrant1(pos,np,ix,iy,x_h,y_h);
  d_right=d2(number(np,iy),number(right_neighbor(np),iy));
  d_upper=d2(number(np,iy),number(upper_neighbor(np),iy));
  f(yrange,xrange) = (iX2(yrange,xrange) ./ d_right) ...
      + (jY2(yrange,xrange) ./ d_upper);
  dfdx(yrange,xrange) = iX(yrange,xrange) ./ d_right;
  dfdy(yrange,xrange) = jY(yrange,xrange) ./ d_upper;
end
if(~(on_left_border(np,ix) || on_upper_border(np,iy)))
  % compute second quadrant
%    fprintf('computing second quadrant...\n')
  [yrange,xrange]=quadrant2(pos,np,ix,iy,x_h,y_h);
  d_left=d2(number(np,iy),number(left_neighbor(np),iy));
  d_upper=d2(number(np,iy),number(upper_neighbor(np),iy));
  f(yrange,xrange) = (iX2(yrange,xrange) ./ d_left) ...
      + (jY2(yrange,xrange) ./ d_upper);
  dfdx(yrange,xrange) = iX(yrange,xrange) ./ d_left;
  dfdy(yrange,xrange) = jY(yrange,xrange) ./ d_upper;
end
if(~(on_left_border(np,ix) || on_lower_border(np,iy)))
  % compute third quadrant
%    fprintf('computing third quadrant...\n')
  [yrange,xrange]=quadrant3(pos,np,ix,iy,x_h,y_h);
  d_left=d2(number(np,iy),number(left_neighbor(np),iy));
  d_lower=d2(number(np,iy),number(lower_neighbor(np),iy));
  f(yrange,xrange) = (iX2(yrange,xrange) ./ d_left) ...
      + (jY2(yrange,xrange) ./ d_lower);
  dfdx(yrange,xrange) = iX(yrange,xrange) ./ d_left;
  dfdy(yrange,xrange) = jY(yrange,xrange) ./ d_lower;
end
if(~(on_right_border(np,ix) || on_lower_border(np,iy)))
  % compute fourth quadrant
%    fprintf('computing fourth quadrant...\n')
  [yrange,xrange]=quadrant4(pos,np,ix,iy,x_h,y_h);
  d_right=d2(number(np,iy),number(right_neighbor(np),iy));
  d_lower=d2(number(np,iy),number(lower_neighbor(np),iy));
  f(yrange,xrange) = (iX2(yrange,xrange) ./ d_right) ...
      + (jY2(yrange,xrange) ./ d_lower);
  dfdx(yrange,xrange) = iX(yrange,xrange) ./ d_right;
  dfdy(yrange,xrange) = jY(yrange,xrange) ./ d_lower;
end

f=hn.*exp(f./(-2));
dfdx=f.*dfdx; dfdy=f.*dfdy;

function b=on_left_border(np,ix)
b = np(1) == ix(1);
function b=on_right_border(np,ix)
b = np(1) == ix(2);
function b=on_lower_border(np,iy)
b = np(2) == iy(1);
function b=on_upper_border(np,iy)
b = np(2) == iy(2);

% compute number of neuron from position
function num=number(np,iy)
ydim=iy(2)-iy(1)+1;
num = np(1)*ydim+abs(np(2)-ydim);

function p=left_neighbor(np)
p = np - [1 0];
function p=right_neighbor(np)
p = np + [1 0];
function p=lower_neighbor(np)
p = np - [0 1];
function p=upper_neighbor(np)
p = np + [0 1];

function [yrange,xrange]=quadrant1(pos,np,ix,iy,x_h,y_h)
if(~on_lower_border(np,iy))
  if(~on_left_border(np,ix))
    yrange=[pos(2)+1:y_h]; xrange=[pos(1)+1:x_h];
  else
    yrange=[pos(2)+1:y_h]; xrange=[1:x_h];
  end
else
  if(~on_left_border(np,ix))
    yrange=[1:y_h]; xrange=[pos(1)+1:x_h];
  else
    yrange=[1:y_h]; xrange=[1:x_h];
  end
end

function [yrange,xrange]=quadrant2(pos,np,ix,iy,x_h,y_h)
if(~on_lower_border(np,iy))
  if(~on_right_border(np,ix))
    yrange=[pos(2)+1:y_h]; xrange=[1:pos(1)];
  else
    yrange=[pos(2)+1:y_h]; xrange=[1:x_h];
  end
else
  if(~on_right_border(np,ix))
    yrange=[1:y_h]; xrange=[1:pos(1)];
  else
    yrange=[1:y_h]; xrange=[1:x_h];
  end
end

function [yrange,xrange]=quadrant3(pos,np,ix,iy,x_h,y_h)
if(~on_upper_border(np,iy))
  if(~on_right_border(np,ix))
    yrange=[1:pos(2)]; xrange=[1:pos(1)];
  else
    yrange=[1:pos(2)]; xrange=[1:x_h];
  end
else
  if(~on_right_border(np,ix))
    yrange=[1:y_h]; xrange=[1:pos(1)];
  else
    yrange=[1:y_h]; xrange=[1:x_h];
  end
end

function [yrange,xrange]=quadrant4(pos,np,ix,iy,x_h,y_h)
if(~on_upper_border(np,iy))
  if(~on_left_border(np,ix))
    yrange=[1:pos(2)]; xrange=[pos(1)+1:x_h];
  else
    yrange=[1:pos(2)]; xrange=[1:x_h];
  end
else
  if(~on_left_border(np,ix))
    yrange=[1:y_h]; xrange=[pos(1)+1:x_h];
  else
    yrange=[1:y_h]; xrange=[1:x_h];
  end
end
