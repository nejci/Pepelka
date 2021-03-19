function [num_cc,cc,cc_idx]=connected_components(adj)
% [num_cc,cc,cc_idx]=connected_components(adj)
%
% connected_components computes the connected components of a graph
% given by adjacency matrix adj. In the case of undirected graphs
% it is expected that adj is symmetric, e.g. adj(k,l)=adj(l,k).
%
% Inputs:
% adj - (matrix) adjacency matrix
% 
% Outputs:
% num_cc - (int) number of connected components in graph
% cc - (cell array of matrices) adjacency matrices of the
%      subgraphs, defined by the connected components
% cc_idx - (cell array of vectors) cc_idx{k} contains the indices
%      of all elements ind adj which are contained in cc{k}, e.g.
%      cc{k}=adj(cc_idx{k},cc_idx{k})
% 
% $Id: connected_components.m 566 2005-04-13 15:50:49Z dome $
% D. Brugger, 12 April 2005
% algo/connected_components.m

if(nargin == 0)
  test_connected_components();
  return;
end

% total number of nodes
n=size(adj,1);
N=1:n;
% number of connected components
num_cc=0;
cc=cell(num_cc,1); cc_idx=cell(num_cc,1);
while(~isempty(N))
  S=reachable_from(adj,N(1));
  N=setdiff(N,S);
  num_cc = num_cc + 1;
  cc{num_cc}=adj(S,S);
  cc_idx{num_cc}=S;
end

function test_connected_components()
% Test case #1 - residual graph of test case #2 in max_flow with
% but regarded as an undirected graph w/r to forward edges
adj=[0 2 9 6 0 0 0; ...
     2 0 0 0 0 0 0; ...
     9 0 0 0 0 0 0; ...
     6 0 0 0 0 0 0; ...
     0 0 0 0 0 22 0; ...
     0 0 0 0 22 0 12; ...
     0 0 0 0 0 12 0];
ecc{1}=[0 2 9 6; ...
        2 0 0 0; ...
        9 0 0 0; ...
        6 0 0 0];
ecc{2}=[0 22 0; ...
        22 0 12; ...
        0 12 0];
ecc_idx{1}=[1 2 3 4];
ecc_idx{2}=[5 6 7];
enum_cc=2;
[num_cc,cc,cc_idx]=connected_components(adj)
check_equal(enum_cc,num_cc,'enum_cc','num_cc');
if(~check_matching(ecc,cc))
  error('check_matching failed\n');
end
if(~check_matching(ecc_idx,cc_idx))
  error('check_matching failed\n');
end

% Test case #2 - residual graph of test case #3 in max_flow 
% but regarded as an undirected graph w/r to forward edges
adj=[0 0 0 0 0 0 0 0 0; ...
     0 0 0 0 0 0 0 0 0; ...
     0 0 0 1 0 0 0 0 0; ...
     0 0 1 0 0 0 0 0 0; ...
     0 0 0 0 0 2 0 1 0; ...
     0 0 0 0 2 0 0 0 0; ...
     0 0 0 0 0 0 0 0 3; ...
     0 0 0 0 1 0 0 0 5; ...
     0 0 0 0 0 0 3 5 0];
enum_cc=4;
ecc{1}=[0]; ecc{2}=[0];
ecc{3}=[0 1; ...
        1 0];
ecc{4}=[0 2 0 1 0; ...
        2 0 0 0 0; ...
        0 0 0 0 3; ...
        1 0 0 0 5; ...
        0 0 3 5 0];
ecc_idx{1}=[1];
ecc_idx{2}=[2];
ecc_idx{3}=[3 4];
ecc_idx{4}=[5 6 7 8 9];
[num_cc,cc,cc_idx]=connected_components(adj)
check_equal(enum_cc,num_cc,'enum_cc','num_cc');
if(~check_matching(ecc,cc))
  error('check_matching failed\n');
end
if(~check_matching(ecc_idx,cc_idx))
  error('check_matching failed\n');
end

fprintf('test_connected_components succeded\n');
