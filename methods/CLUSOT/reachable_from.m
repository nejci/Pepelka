function [S]=reachable_from(adj,s)
% [S]=reachable_from(adj,s)
%
% reachable_from determines all nodes in the graph given by the
% adjacency matrix adj which are reachable from node s
%
% Inputs:
% adj - (matrix) adjacency matrix of graph
% s - (int) start node index
%
% Outputs:
% S - (vector) set of all nodes reachable from node s
% 
% $Id: reachable_from.m 564 2005-04-13 15:46:12Z dome $
% D. Brugger, 12 April 2005
% algo/reachable_from.m

if(nargin == 0)
  test_reachable_from();
  return;
end

S=reachable_from_helper(adj,s,[]);

% determine which nodes are reachable from s, e.g. compute
% the transitive closure of s w/r to the graph given by adj
function S=reachable_from_helper(adj,s,S)
if(isempty(S))
  S=[s];
end
n=size(adj,2);
mask=ones(1,n); mask(S)=0;
ind=find((adj(s,:).*mask) > 0);
%ind=setdiff(find(adj(s,:) > 0),S);
if(~isempty(ind))
  S=union(S,ind);
  for k=1:size(ind,2)
    S=union(S,reachable_from_helper(adj,ind(k),S));
  end
end

function test_reachable_from()
% Test case #1, residual graph of test case #3 in max_flow_test
adj=[0 0 0 0 0 0 0 0 0; ...
     2 0 0 0 0 0 0 0 0; ...
     2 0 0 1 0 0 0 0 0; ...
     0 2 3 0 0 0 0 0 0; ...
     0 0 0 4 0 2 0 1 0; ...
     0 0 0 0 4 0 0 0 0; ...
     0 0 0 0 2 0 0 0 3; ...
     0 0 0 0 1 0 0 0 5; ...
     0 0 0 0 0 2 5 5 0];
eS=[1];
S=reachable_from(adj,1)
check_equal(eS,S,'eS','S');
% Test case #2, residual graph of test case #2 in max_flow_test
adj=[0 2 9 6 0 0 0; ...
    18 0 0 0 0 0 0; ...
     6 5 0 0 0 0 0; ...
     4 4 0 0 0 0 0; ...
     0 5 0 8 0 22 0; ...
     0 9 6 0 8 0 12; ...
     0 0 0 0 10 18 0];
eS=[1,2,3,4];
S=reachable_from(adj,1)
check_equal(eS,S,'eS','S');

fprintf('test_reachable_from succeded\n');