function D = sqdistance2(A,B)
% Compute squared Euclidean distances between all pair of points.
%   A: N1 x D data matrix
%   B: N2 x D data matrix
%   D: N1 x N2 pairwise square distance matrix
% Written by Michael Chen (sth4nth@gmail.com).
% Modified by Nejc Ilc

if nargin == 1
    A = bsxfun(@minus,A,mean(A,1));
    S = full(dot(A,A,2));
    D = bsxfun(@plus,S,S')-full(2*(A*A'));
elseif nargin == 2
    assert(size(A,2)==size(B,2));
    D = bsxfun(@plus,full(dot(B,B,2))',full(dot(A,A,2)))-full(2*(A*B'));
end

D = max(0,D); % avoid accuracy errors with negative distance
