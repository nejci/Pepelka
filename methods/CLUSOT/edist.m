function [d]=edist(v1,v2)
% [d]=edist(v1,v2)
%
% edist computes the euclidian distance between vectors v1 and v2.
%
% Inputs:
% v1,v2 - (vector)
% 
% Output:
% d - (double)
%
% $Id: edist.m 488 2005-03-25 12:43:03Z dome $
% D. Brugger, 25 March 2005
% util/edist.m

if(nargin == 0)
  test_edist();
  return;
end

d=sqrt(sum((v1-v2).^2));

function test_edist()
v1=[0 1]; v2=[1 1];
ed=1;
d=edist(v1,v2)
check_equal(d,ed,'d','ed');

v1=[1 2]; v2=[-2 2];
ed=3;
d=edist(v1,v2)
check_equal(d,ed,'d','ed');

v2=[-2 3];
ed=sqrt(10);
d=edist(v1,v2)
check_equal(d,ed,'d','ed');

fprintf('test_edist succeded\n');