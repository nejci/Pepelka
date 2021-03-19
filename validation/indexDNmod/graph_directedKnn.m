function [G,W] = graph_directedKnn(M,k,which_matrix)
% GRAPH_DIRECTEDKNN creates directed kNN graph on M
% Usage: [G,W] = graph_directedKnn(M,k,which_matrix) 
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
% W             : adjacency matrix of the directed knn graph
%--------------------------------------------------------------------------
% For a similarity matrix M, returns the directed knn graph, edges are
% weighted by S. For a distance matrix D, returns the directed knn graph,
% edges are weighted by distance. Self-edges are excluded in both graphs.
%--------------------------------------------------------------------------
% Copyright 2007, Matthias Hein and Ulrike von Luxburg.
%--------------------------------------------------------------------------

% implemented by brute force sorting
%
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

% testing whether matrix is square: 
if (size(M,1) ~= size(M,2))
  error('Matrix not square!')
end


n = size(M,1); 

% to exclude self-edges, set diagonal of sim/dissim matrix to Inf or 0
if (strcmpi(which_matrix,'sim'))
    M(1:n+1:n^2) = 0;
elseif (strcmpi(which_matrix, 'dist'))
    M(1:n+1:n^2) = Inf;
else
    error('Unknown matrix type')
end



% now do it: 
W = M; 

if (strcmp(which_matrix, 'sim'))
  
  for it = 1:n
    % sort points according to similarities: 
    [~,order] = sort(M(it,:), 'descend'); 
    
    % for all points which are not among the k nearest neighbors, set W to 0: 
    W(it, order(k+1:end)) = 0;
    %W(it,order(1:k) just stays the same
  end
  
elseif (strcmp(which_matrix, 'dist'))
  D=M;
  for it = 1:n
    % sort points according to distances: 
    [~,order] = sort(M(it,:), 'ascend'); 
    % for all points which are not among the k nearest neighbors, set W to 0: 
    W(it, order(k+1:end)) = 0;
    W(it, order(1:k)) = 1; %unweighted!
	D(it, order(k+1:end)) = 0;
  end
  
else 
  error('build_directed_knn_graph: unknown matrix type')
end

% default MATLAB format for storing graphs is sparse matrix
G=sparse(D);


  
