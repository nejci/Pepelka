function [r]=bresenham(a,b)
% [r]=bresenham(a,b)
%
% bresenham computes a rasterized line (set of points with integer 
% coordinates) for a line through points a and b.
%
% Inputs:
% a, b - (vector) line endpoints with integer values, e.g. in Z^2
%
% Output:
% r - (matrix) raster points along line given by a and b,
%      a n x 2 matrix, where n is the number of raster points.
% 
% $Id: bresenham.m 408 2005-03-03 11:24:51Z dome $
% D. Brugger, 27 February 2005
% algo/bresenham.m

if(nargin == 0)
  test_bresenham();
  return;
end

a_x=a(1); a_y=a(2); b_x=b(1); b_y=b(2); 
dx = b_x - a_x; dy = b_y - a_y;
adx = abs(dx); ady = abs(dy);

% handle special cases dx=0 or dy=0
if(dx == 0)
  if(dy < 0)
    inc_y = -1;
  else
    inc_y = 1;
  end
  r=zeros(ady+1,2);
  for k=1:ady+1
    r(k,:)=[a_x a_y];
    a_y = a_y + inc_y;
  end
  return;
end
if(dy == 0)
  if(dx < 0)
    inc_x = -1;
  else
    inc_x = 1;
  end
  r=zeros(adx+1,2);
  for k=1:adx+1
    r(k,:)=[a_x a_y];
    a_x = a_x + inc_x;
  end
  return;
end

swap_flag=0; m=dy/adx;
% swap x and y if we are in II,III,VI or VII octant
if(m > 1 || m < -1)
  [a_x,a_y]=swap(a_x,a_y); [b_x,b_y]=swap(b_x,b_y);
  swap_flag=1;
  % recalculate dx,dy,adx,ady
  dx = b_x - a_x; dy = b_y - a_y;
  adx = abs(dx); ady = abs(dy);
end

% check if we are in I,IV,V,VIII octant and set 
% increments accordingly
if(dx > 0)
  inc_x = 1;
else
  inc_x = -1;
end
if(dy > 0)
  inc_y = 1;
else
  inc_y = -1;
end

x=a_x; y=a_y;
r=zeros(adx+1,2); 
% insert starting point
r(1,:)=[a_x a_y]; x = x + inc_x;
p=2*ady - adx; dp1=2*(ady-adx); dp2=2*ady;
for k=2:adx+1
  if(p > 0) % next point is right and up
    y = y + inc_y;
    r(k,:)=[x,y];
    p = p + dp1;
  else % next point is right
    r(k,:)=[x,y];
    p = p + dp2;
  end
  x = x + inc_x;
end

if(swap_flag)
  r = flipdim(r,2);
end

function [x,y]=swap(x,y)
tmp = x; x = y; y = tmp;

function test_bresenham()
% Test case #0: long line
a=[4 11]; b=[19 8];
er=[4 11; ...
    5 11; ...
    6 11; ...
    7 10; ...
    8 10; ...
    9 10; ...
    10 10; ...
    11 10; ...
    12 9; ...
    13 9; ...
    14 9; ...
    15 9; ...
    16 9; ...
    17 8; ...
    18 8; ...
    19 8];
r=bresenham(a,b)
check_equal(er,r,'er','r');

% Test case #1: I & II octant
a=[0 0]; b1=[5 2]; b2=[1 4];
er1=[0 0; ...
     1 0; ...
     2 1; ...
     3 1; ...
     4 2; ...
     5 2];
er2=[0 0; ...
     0 1; ...
     0 2; ...
     1 3; ...
     1 4];
r1=bresenham(a,b1)
check_equal(er1,r1,'er1','r1');
r2=bresenham(a,b2)
check_equal(er2,r2,'er2','r2');

% Test case #2: III & IV octant
a1=[-1 0]; b1=[-5 3]; 
a2=[0 2]; b2=[-2 6];
er1=[-1 0; ...
     -2 1; ...
     -3 1; ...
     -4 2; ...
     -5 3];
er2=[0 2; ...
     0 3; ...
     -1 4; ...
     -1 5; ...
     -2 6];
r1=bresenham(a1,b1)
check_equal(er1,r1,'er1','r1');
r2=bresenham(a2,b2)
check_equal(er2,r2,'er2','r2');

% Test case #3: V & VI octant
a=[0 0]; b1=[-5 -2]; b2=[-2 -6];
er1=[0 0; ...
     -1 0; ...
     -2 -1; ...
     -3 -1; ...
     -4 -2; ...
     -5 -2];
er2=[0 0; ...
     0 -1; ...
     -1 -2; ...
     -1 -3; ...
     -1 -4; ...
     -2 -5; ...
     -2 -6];
r1=bresenham(a,b1)
check_equal(er1,r1,'er1','r1');
r2=bresenham(a,b2)
check_equal(er2,r2,'er2','r2');

% Test case #4: VII & VIII octant
a=[0 0]; b1=[6 -3]; b2=[3 -4];
er1=[0 0; ...
     1 0; ...
     2 -1; ...
     3 -1; ...
     4 -2; ...
     5 -2; ...
     6 -3];
er2=[0 0; ...
     1 -1; ...
     1 -2; ...
     2 -3; ...
     3 -4];
r1=bresenham(a,b1)
check_equal(er1,r1,'er1','r1');
r2=bresenham(a,b2)
check_equal(er2,r2,'er2','r2');

% Test case #5: dx=0 | dy=0
a=[0,0]; b1=[0 5]; b2=[2 0]; b3=[0 -3]; b4=[-2 0];
er1=[0 0; ...
     0 1; ...
     0 2; ...
     0 3; ...
     0 4; ...
     0 5];
er2=[0 0; ...
     1 0; ...
     2 0];
er3=[0 0;
     0 -1;
     0 -2;
     0 -3];
er4=[0 0; ...
     -1 0; ...
     -2 0];
r1=bresenham(a,b1)
check_equal(er1,r1,'er1','r1');
r2=bresenham(a,b2)
check_equal(er2,r2,'er2','r2');
r3=bresenham(a,b3)
check_equal(er3,r3,'er3','r3');
r4=bresenham(a,b4)
check_equal(er4,r4,'er4','r4');

fprintf('test_bresenham succeded!\n')

% put this into test/check_equal.m
%function check_equal(er, r)
%[n1,m1]=size(er); [n2,m2]=size(r);
%if(n1 ~= n2 || m1 ~= m2)
%  error(['!!! Result size not equal: size(er)=[%d,%d] != size(r)=[%d,%d] ' ...
%         ' !!!'], n1, m1, n2, m2);
%end
%if(sum(sum(er == r)) ~= n1*m1)
%  error('!!! Result r does not match expected result er !!!');
%end
