function [GC, TR] = graph_globalCC(G)
%
% Computes global clustering coefficient and transitivity ratio of graph G.
% G is a symmetric square adjancency matrix representing an undirected graph
%
% GC = 3 * number of triangles / number of all triples
% TR = 3 * number of triangles / number of connected triples
%
% Written by Nejc Ilc, 2013.
% Using code from: 
%   Matlab Tools for Network Analysis
%   Copyright (c) 2011, Massachusetts Institute of Technology.

if ~isequal(G,G'), error('Only undirected (symmetric) graphs allowed.'); end

G = G > 0;  % binarize, no multiple edges

% Compute the number of triangles (3-loops)
numTriangles = trace(G^3)/6;

% Compute the number of all triples
numNeigh = sum(G,2);
numC = pplk_nchoosek2(numNeigh,2);
numC(isnan(numC)) = 0;
triplesAll = sum(numC);
% SLOWER ALTERNATIVE
% triplesAll = 0;  % initialize
% for i=1:length(G)
%     neigh = find(G(i,:)>0);
%     numNeigh = length(neigh);
%     % handle leaves, no triple here
%     if numNeigh < 2;
%         continue;
%     end
%     triplesAll = triplesAll + nchoosek(numNeigh,2);
% end


GC = 3*numTriangles / triplesAll;

if nargout > 1
    % due to the symmetry triangles repeat 3 times in the nchoosek count
    triplesConn = triplesAll - 2*numTriangles;
    TR = 3*numTriangles / triplesConn;
end

