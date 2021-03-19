function [GCw] = graph_globalCCw(G,method,isDist,dist2simAlg)
% Weighted global clustering coefficient 
% as proposed by Opsahl & Panzarasa, 2009
%--------------------------------------------------------------------------
% G:      symmetric weight matrix of a graph (undirected)
% method: how to compute the triplet value:
%         ['am'] - arithmetic mean
%          'gm'  - geometric mean
%          'ma'  - maximum weight
%          'mi'  - minimum weight
%
% isDist: [0] - weigths are similarities (do not transform) 
%          1  - weights are distances (needs transformation to similarities),
%
% dist2simAlg: algorithmm to transform distances to similarities:
%              'lin'    - linear, normalized
%              ['lin2'] - linear, normalized, max=1, min sim. = min. dist
%              'gauss'  - gaussian; sigma is half of range if not otherwise defined!
%              'gauss2' - quasi gaussian, normalized, without sigma
%              'gauss3' - quasi gaussian, without sigma
%--------------------------------------------------------------------------
% References: 
% T. Opsahl and P. Panzarasa, �Clustering in weighted networks,� Social
%   Networks, vol. 31, no. 2, pp. 155�163, May 2009.
%--------------------------------------------------------------------------
% Author: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
% Last modified: 25-November-2013 by Nejc Ilc
%--------------------------------------------------------------------------

n = length(G);
if ~isequal(G,G'), error('Only (undirected) symmetric graphs allowed'); end

if ~exist('method','var') || isempty(method)
    method = 'am';
end
switch method
    case 'am'
        fun = @(Pw) mean(Pw,2);
    case 'gm'
        fun = @(Pw) sqrt(prod(Pw,2));
    case 'mi'
        fun = @(Pw) min(Pw,[],2);
    case 'ma'
        fun = @(Pw) max(Pw,[],2);
    otherwise
        error('Wrong method! Select from: ''am'', ''gm'', ''mi'', ''ma''.');
end

if ~exist('isDist','var') || isempty(isDist)
    isDist = 0;
end
if ~exist('dist2simAlg','var') || isempty(dist2simAlg)
    dist2simAlg = 'lin2';
end
% transform distances to similarities, which are used as weight in a graph
if isDist
    G = dist2sim(G,dist2simAlg);
end

% delete self-loops
G = G .* ~eye(n);

sumTriangles = 0;
sumTriples = 0;

for i = 1:n
    neigh = find(G(:,i));
    a = length(neigh); % degree of a node i
    if a < 2
        continue;
    end

    neighW = G(i,neigh);
    neighPw = VChooseK(neighW,2);
    w = fun(neighPw);
    sumTriples = sumTriples + sum(w);
    triInd = tril(G(neigh,neigh)) > 0;
    indMat = bsxfun(@minus,reshape(1:a^2,a,a),cumsum(1:a));
    isTriangle = indMat(triInd);
    sumTriangles = sumTriangles + sum(w(isTriangle));
    
    % SLOWER ALTERNATIVE
    %     for k1 = 1:(a-1)
    %         for k2 = k1+1:a
    %             w = sqrt( G(i,neigh(k1)) * G(i,neigh(k2)) );
    %             % consider triplet
    %             sumTriples = sumTriples + w;
    %             % is there a triangle?
    %             if G(neigh(k1),neigh(k2))
    %                 sumTriangles = sumTriangles + w;
    %             end
    %         end
    %     end
    
end
GCw = sumTriangles/sumTriples;
