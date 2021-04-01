function [G,W] = graph_symmetricKnn(M,k,which_matrix)
% Usage:  W = GD_BuildSymmetricKnnGraph(M,k,which_matrix)
%
% Input: 
% M                = either the distance or the similarity matrix, needs to be square
% k                = parameter of the knn graph
% which_matrix     = either 'sim' or 'dist'
% 
% Output: 
% G				 = [n x n] sparse matrix that represents a graph. Nonzero
%					entries in matrix G represent the weights of the edges.
% W              = adjacency matrix of the symmetric knn graph without self-edges                 
%
% For a similarity matrix S, returns the undirected knn graph, edges are weighted by S. 
% For a distance matrix D, returns the undirected (unweighted!) knn graph. If you want to get 
% a weighted graph in this case, you need to take care of transforming D to S yourself and then 
% call the function with a similarity matrix 
%
% Excludes self-edges. 

%implemented brute force 



if (size(M,1) ~= size(M,2))
  error('Matrix not square!')
end

% compute directed KNN graph: 
[dummy,W] = graph_directedKnn(M,k,which_matrix); 
% transform it to symmetric one: 
W = max(W,W'); 

% weight by distances
if strcmp(which_matrix,'dist')
	G=W.*M;
	G=sparse(tril(G));
end
