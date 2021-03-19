function [G,W] = graph_mutualKnn(M,k,which_matrix)
% GRAPH_MUTUALKNN creates mutual kNN graph on M
% Usage:  [G,W] = graph_mutualKnn(M,k,which_matrix)
%--------------------------------------------------------------------------
% INPUTS: 
% M             : either the distance or the similarity matrix, needs to be 
%                 square, symmetric, non-negative
% k             : connectivity parameter of the kNN graph
% which_matrix  : either 'sim' or 'dist' (similarity or distance matrix)
%--------------------------------------------------------------------------
% OUTPUTS: 
% G				: [n x n] sparse matrix that represents a graph. Nonzero
%				  entries in matrix G represent the weights of the edges.
% W             : adjacency matrix of the mutual kNN graph
%--------------------------------------------------------------------------
%
% For a similarity matrix S, returns the mutual knn graph, edges are
% weighted by S. For a distance matrix D, returns the undirected
% (unweighted!) mutual  knn graph. If you want to get a weighted graph in
% this case, you need to take care of transforming D to S yourself and then
% call the function with a similarity matrix.
% Self edges are excluded in both cases. 
%--------------------------------------------------------------------------
% Copyright 2007, Matthias Hein and Ulrike von Luxburg.
%--------------------------------------------------------------------------

% Written by Matthias Hein and Ulrike von Luxburg. 
% It has been first used at the practical sessions of the Machine Learning 
% Summer School 2007 at the Max Planck Institute for Biological Cybernetics 
% in Tuebingen, Germany: http://www.mlss.cc/tuebingen07/
% 
% This file is released as free software. It can be used by everybody, as
% long as this file and the credit comments in the demo are not removed.
% 
% Please send comments, suggestions and bug reports to:
% Matthias Hein
% Saarland University, Germany
% http://www.ml.uni-saarland.de/contact.html
%            or 
% Ulrike von Luxburg
% Max Planck Institute for Biological Cybernetics, Tuebingen, Germany
% http://www.kyb.mpg.de/~ule

% check: 
if (size(M,1) ~= size(M,2))
  error('Matrix not square!')
end

% build the directed knn graph: 
[dummy,W] = graph_directedKnn(M,k,which_matrix); 

% transform it to the mutual one: 
W = min(W,W');

% weight by distances
if strcmp(which_matrix,'dist')
	G=W.*M;
	G=sparse(tril(G));
end