function [CCw,avgCCw]=graph_localCCw(G,method,isDist,dist2simAlg,nanMode)

% Weighted local clustering coefficient
%--------------------------------------------------------------------------
% G: weight matrix of a graph (directed or undirected)
%
% method: considered only in case of undirected graph.
%         ['onnela'] - using method by Onnela et al., 2005
%          'zhang'   - using method by Zhang & Horvath, 2005
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
%
% nanMode:   
%           [0] - the CC of node with degree <= 1 is 0 and is included in
%                 the average.                 
%            1  - the CC of node with degree <= 1 is NaN and is not
%                 included in the average.
%
%--------------------------------------------------------------------------
% Weighting algorithm for undirected graphs adopted by Onnela et al. 2005.
% Weighting algorithm for directed graphs adopted by Fagiolo 2007.
% Code inspired by Mikail Rubinov (Brain Connectivity Toolbox)
%--------------------------------------------------------------------------
% References:
% Complex network measures of brain connectivity: Uses and interpretations.
%   Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.
%
% J.-P. Onnela, J. Saram�ki, J. Kert�sz, and K. Kaski, �Intensity and
%   coherence of motifs in weighted complex networks,� Physical Review E,
%   vol. 71, no. 6, p. 065103, Jun. 2005.
%
% B. Zhang and S. Horvath, �A general framework for weighted gene
%   co-expression network analysis,� Statistical applications in genetics
%   and molecular biology, vol. 4, no. 1, p. Article17, Jan. 2005.
%
% G. Fagiolo, �Clustering in complex directed networks,� Physical Review E,
%   vol. 76, no. 2, p. 026107, Aug. 2007.
%--------------------------------------------------------------------------
% Author: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
% Last modified: 25-November-2013 by Nejc Ilc
%--------------------------------------------------------------------------

n = length(G);
if isequal(G,G')
    isdirected = 0;
else
    isdirected = 1;
end
if ~isdirected && (~exist('method','var') || isempty(method))
    method = 'onnela';
end
if ~exist('isDist','var') || isempty(isDist)
    isDist = 0;
end
if ~exist('dist2simAlg','var') || isempty(dist2simAlg)
    dist2simAlg = 'lin2';
end
if ~exist('nanMode','var') || isempty(nanMode)
    nanMode = 0;
end

% transform distances to similarities, which are used as weight in a graph
if isDist
    % linear
    G = dist2sim(G,dist2simAlg);
end

% delete self-loops
G = G .* ~eye(n);
% normalize weigths by the maximum weight
G = G./max(G(:));     

if isdirected
    A = W~=0;                     % adjacency (binary) matrix
    S = W.^(1/3)+(W.').^(1/3);	  % symmetrized weights matrix ^1/3
    K = sum(A+A.',2);             % total degree (in + out)
    cyc3 = diag(S^3)/2;           % number of 3-cycles (ie. directed triangles)
    K(cyc3==0) = Inf;             % if no 3-cycles exist, make C=0 (via K=inf)
    CYC3 = K.*(K-1)-2*diag(A^2);  % number of all possible 3-cycles
    CCw = cyc3./CYC3;             % clustering coefficient
else
    if strcmp(method,'onnela')
        % Method by Onnela et al.
        K = sum(G~=0,2);          % degrees of nodes
        nanMask = false(n,1);
        if nanMode
            nanMask = K<2;
        end
        cyc3 = diag((G.^(1/3))^3);  % values of triangles
        K(cyc3==0) = Inf;           % if no triangles exist, make C=0 (via K=inf)
        CCw = cyc3./(K.*(K-1));     % clustering coefficient
        CCw(nanMask) = NaN;
        
    elseif strcmp(method,'zhang')
        % Method by Zhang & Horvath
        % Implemented as pointed out in Eq. 4 in paper:
        % G. Kalna and D. J. Higham, �A clustering coefficient for weighted
        % networks, with application to gene expression data,� AI Communications,
        % vol. 20, no. 4, pp. 263�271, Jan. 2007.
        W3kk = diag(G^3);      % weighted value of triangles        
        sumSqWk = sum(G.^2,1); % sum of squared weights
        sqSumWk = sum(G,1).^2; % squared sum of weights        
        CCw = W3kk ./ (sqSumWk - sumSqWk)';
        if ~nanMode
            CCw(isnan(CCw)) = 0;
        end
    else
        error('Wrong method!');
    end   
end

% SLOW alternative for Onnela's method
% CCw = zeros(1,n);
% for i = 1:n
%     neighbours = find(G(:,i))';
%     a = length(neighbours); % degree of a node i
%     if a < 2
%         continue;
%     end
%     E = 0;
%     for k1 = 1:(a-1)
%         for k2 = k1+1:a
%             if G(neighbours(k1),neighbours(k2)) %|| A(neighbours(k2),neighbours(k1))
%                 % obtain weights in the triangle (k,k1,k2)
%                 w = ( G(i,neighbours(k1)) * G(neighbours(k1),neighbours(k2)) * G(i,neighbours(k2)) )^(1/3);
%                 E = E + w;
%             end
%         end
%     end
%     CCw(i) = 2 * E / (a * (a-1));
% end

if nargout > 1
    CCv = CCw(~isnan(CCw));
    avgCCw = mean(CCv);
end
