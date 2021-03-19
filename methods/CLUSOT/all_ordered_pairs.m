function [pairs]=all_ordered_pairs(v)
% [pairs]=all_ordered_pairs(v)
%
% all_ordered_pairs returns a matrix of all ordered pairs of
% elements of vector v
%
% Input:
% v - (vector)
%
% Output:
% pairs - (matrix)
% 
% $Id: all_ordered_pairs.m 566 2005-04-13 15:50:49Z dome $
% D. Brugger, 13 April 2005
% algo/all_ordered_pairs.m

if(nargin == 0)
  test_all_ordered_pairs();
  return;
end

% handle special case
if(size(v,2) == 1)
  pairs=[];
end

pos = 1;
for k=1:size(v,2)
  for l=k+1:size(v,2)
    pairs(pos,:)=[k l];
    pos = pos + 1;
  end
end

function test_all_ordered_pairs()
v = [1:5];
epairs=[1 2; ...
        1 3; ...
        1 4; ...
        1 5; ...
        2 3; ...
        2 4; ...
        2 5; ...
        3 4; ...
        3 5; ...
        4 5];
pairs=all_ordered_pairs(v)
check_equal(epairs,pairs,'epairs','pairs');
v=[1];
epairs=[];
pairs=all_ordered_pairs(v)
check_equal(epairs,pairs,'epairs','pairs');

fprintf('test_all_ordered_pairs succeded\n');
