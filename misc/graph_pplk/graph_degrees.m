function [deg,indeg,outdeg]=graph_degrees(G,isdirected)
% Compute the total degree, in-degree and out-degree of a graph based on
% the adjacency matrix; should produce weighted degrees, if the input matrix is weighted
% INPUTS: adjacency matrix
% OUTPUTS: degree, indegree and outdegree sequences
% GB, Last Updated: October 2, 2009
% Mod by Nejc Ilc

indeg = sum(G,1);
outdeg = sum(G,2);

if isdirected
    % total degree
    deg = indeg + outdeg;
else
    % undirected graph: indeg=outdeg
    % add self-loops twice, if any
    deg = indeg + diag(G)';
end