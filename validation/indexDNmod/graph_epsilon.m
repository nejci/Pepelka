function [G,W] = graph_epsilon(M,t,which_matrix)
% Usage: [G,W] = graph_epsilon(M,t,which_matrix)
%
% Input:
% M                = either distance matrix or similarity matrix, needs to be square and symmetric
% t                = threshold of the eps graph
% which_matrix     = 'sim' or 'dist'
%
% Output:
% G				 = [n x n] sparse matrix that represents a graph. Nonzero
%					entries in matrix G represent the weights of the edges.
% W              = adjacency matrix of the epsilon-neighborhood (UNDIRECTED) graph
%
% For a distance matrix: connects all points with distanace M(it,jt) <= t.
% For a similarity matrix: connects all points with similarity M(it,jt) >= t, and weights the edges by the similarity.
% Implemented brute force
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

if (size(M,1) ~= size(M,2))
    error('Matrix not square!')
end
n = size(M,1);

if isempty(t)
    t=defaultEps(M);
end

if (strcmpi(which_matrix,'sim'))
    % to exclude self-edges, set diagonal to 0
    M(1:n+1:n^2) = 0;
    W = (M >= t) .* M;
    %  W = +W; % make it numeric rather than logical
    
elseif (strcmpi(which_matrix, 'dist'))
    % to exclude self-edges, set diagonal of sim/dissim matrix to Inf or 0
    M(1:n+1:n^2) = 0;
    W = (M <= t) .* M;
    %W = +W; % make it numeric rather than logical
    
else
    error('Unknown matrix type')
end

% default MATLAB format for storing graphs is sparse matrix
G=sparse(tril(W));


function goodPara = defaultEps(D)
% set current value of current eps/sigma to something useful in the beginning:


if (isempty(D))
    goodPara = 0.5; %default
    
else
    % try to select reasonable eps according to bold rule of thumb
    % mean of similarity of k-th nearest neighbor:
    k = min(7,size(D,1)-1);
    Dsorted = sort(D,2,'ascend');
    Dk = Dsorted(:,k+1); %attention, start to count at second col, first is always 1
    goodPara = mean(Dk);
end
