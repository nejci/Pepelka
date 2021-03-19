function [r]=lin_scale(val,i1,i2)
% [r]=lin_scale(val,i1,i2)
%
% lin_scale linearly scales all values in matrix val
% from interval i1 to interval i2
%
% Inputs:
% val - (matrix) values to be scaled
% i1,i2 - (vectors) old and new interval for values, e.g.
%          i1=[old_min,old_max] and i2=[new_min,new_max]
%
% Outputs:
% r - (matrix) scaled values
% 
% $Id: lin_scale.m 440 2005-03-16 11:57:58Z dome $
% D. Brugger, 16 March 2005
% util/lin_scale.m

if(nargin == 0)
  test_lin_scale();
  return;
end

old_min=i1(1); old_max=i1(2);
new_min=i2(1); new_max=i2(2);
% scaling factor
scf=(new_max-new_min)/(old_max-old_min);
r=(val-old_min).*scf+new_min;

function test_lin_scale()
val=[4 3 2 5];
i1=[2 5]; i2=[2 4];
er=[2+4/3 2+2/3  2 4];
r=lin_scale(val,i1,i2)
check_equal(r,er,'r','er');

val=[1 2 4 2; ...
     8 6 5 3];
er=[2 2+8/7 2+3*8/7  2+8/7
    10 2+5*8/7 2+4*8/7 2+2*8/7]
i1=[1 8]; i2=[2 10];
r=lin_scale(val,i1,i2)
check_equal(r,er,'r','er',0.000001);

fprintf('test_lin_scale succeded\n');
