function [G,W] = graph_epsilon(M,t,which_matrix)
% Usage: W = GD_BuildEpsilonGraph(M,t,which_matrix)
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

%implemented brute force

if (size(M,1) ~= size(M,2))
    error('Matrix not square!')
end
n = size(M,1);

if isempty(t)
    t=defaultEps(M);
end

if (strcmp(which_matrix,'sim'))
    % to exclude self-edges, set diagonal to 0
    for it=1:n
        M(it,it) = 0;
    end    
    W = (M >= t) .* M;    
    %  W = +W; % make it numeric rather than logical
    
elseif (strcmp(which_matrix, 'dist'))
    % to exclude self-edges, set diagonal of sim/dissim matrix to Inf or 0
    for it=1:n
        M(it,it) = 0;
    end
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
    k = 7;
    Dsorted = sort(D,2,'ascend');
    Dk = Dsorted(:,k+1); %attention, start to count at second col, first is always 1
    goodPara = mean(Dk);
end
