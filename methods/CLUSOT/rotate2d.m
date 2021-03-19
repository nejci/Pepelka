function [r]=rotate2d(v,delta)
% [r]=rotate2d(v,delta)
%
% rotate2d rotates vector v counterclockwise by angle delta, given
% in radians, e.g. rotate2d([1 0],pi/2) => [0 1]
%
% Inputs:
% v - (vector)
% delta - (angle)
%
% Ouput:
% r - (vector)
% 
% $Id: rotate2d.m 488 2005-03-25 12:43:03Z dome $
% D. Brugger, 25 March 2005
% util/rotate2d.m

if(nargin == 0)
  test_rotate2d();
  return;
end

% rotation matrix
M = [cos(delta) -sin(delta); ...
     sin(delta) cos(delta)];
if(size(v,1) == 1)
  r = M*v'; r=r';
else
  r = M*v;
end

function test_rotate2d()
v=[1 0]; delta=pi/2;
er=[0 1];
r=rotate2d(v,delta)
check_equal(r,er,'r','er',0.000001);
v=v'; er=er';
r=rotate2d(v,delta)
check_equal(r,er,'r','er',0.000001);

v=[0 -1]; delta=pi;
er=[0 1];
r=rotate2d(v,delta)
check_equal(r,er,'r','er',0.000001);

fprintf('test_rotate2d succeded\n')
