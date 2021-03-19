function [s1,s2]=linseg(M,alpha,ix,iy)
% [s1,s2]=linseg(M,alpha,ix,iy)
%
% linseg computes the two intersection points s1,s2 of linesegments running
% through point M, having a slope given by angle alpha, with the
% borders of a rectangle given by intervals ix=[x_low,x_high] and
% iy=[y_low, y_high]. 
%
% Inputs:
% M - point in Z^2 
% alpha - angle in [0,pi]
% ix - x-interval of bounding rectangle
% iy - y-interval of bounding rectangle
%
% Outputs:
% s1,s2 - intersection points in Z^2
% 
% $Id: linseg.m 409 2005-03-04 19:12:19Z dome $
% D. Brugger, 28 February 2005
% algo/linseg.m

if(nargin == 0)
  test_linseg();
  return;
end

m_x=M(1); m_y=M(2);
x_l=ix(1); x_h=ix(2);
y_l=iy(1); y_h=iy(2);

% deal with special cases
if(alpha == 0 || alpha == pi)
  s1 = [x_l, m_y]; s2 = [x_h, m_y];
  return;
end
if(alpha == pi/2)
  s1 = [m_x, y_l]; s2 = [m_x, y_h];
end

% now we have alpha > 0 and alpha != pi/2 and alpha < pi
m = sin(alpha)/cos(alpha); mi = 1/m;
s1 = []; s2 = []; s3 = []; s4 = [];

x = -mi*(m_y - y_l)+m_x;
if(x >= x_l && x <= x_h)
  s1 = [x, y_l];
end
x = -mi*(m_y - y_h)+m_x;
if(x >= x_l && x <= x_h)
  s2 = [x, y_h];
end
y = m*(x_l - m_x)+m_y;
if(y >= y_l && y <= y_h)
  if(y <= m_y)
    s1 = [x_l, y];
  else
    s2 = [x_l, y];
  end
end
y = m*(x_h - m_x)+m_y;
if(y >= y_l && y <= y_h)
  if(y > m_y)
    s2 = [x_h, y];
  else
    s1 = [x_h, y];
  end
end
% make s1 and s2 have integer coordinates
% by rounding to nearest integer
s1=round(s1); s2=round(s2);

function test_linseg()
m=[2,1]; ix=[0,4]; iy=[0,3];
% test case #1
alpha=0;
es1=[0,1]; es2=[4,1];
[s1,s2]=linseg(m,alpha,ix,iy)
check_equal(s1,s2,es1,es2);
% test case #2
alpha=pi/2;
es1=[2,0]; es2=[2,3];
[s1,s2]=linseg(m,alpha,ix,iy)
check_equal(s1,s2,es1,es2);
% test case #3
alpha=pi/4;
es1=[1,0]; es2=[4,3];
[s1,s2]=linseg(m,alpha,ix,iy)
check_equal(s1,s2,es1,es2);

fprintf('test_linseg succeded!\n')

function check_equal(s1,s2,es1,es2)
if(sum(es1 == s1) ~= 2 || sum(es2 == s2) ~= 2)
  error(['!!! Calculated result s1=(%d,%d), s2=(%d,%d) does not match ' ...
         'expected result es1=(%d,%d), es2=(%d,%d)!!!'], s1(1), ...
        s1(2), s2(1), s2(2), es1(1), es1(2), es2(1), es2(2));
end
