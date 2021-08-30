function S = pplk_dist2sim(D,method,isGraph,sigma,order)
% S = pplk_dist2sim(D,method,isGraph,sigma,order)
% Converts distance matrix D to similarity matrix S, using transformation.
%
% INPUTS
% method 
%   'lin'
%       Linear, normalized. 
%   'lin2'
%       Linear, normalized, max=1, min sim. = min. dist.
%   'ncut'
%       As defined by authors of Normalized Cuts algorithm. 
%       S = exp(-(D./sigma).^order). 
%   'gauss'
%       Gaussian; sigma is half of range if not otherwise defined! 
%   'gauss2'
%       Quasi gaussian, normalized, without sigma.
%   'gauss3'
%       Quasi gaussian, without sigma.
%
% isGraph
%   If D is distance matrix of graph, set this to 1 to ensure zero entries
%   remain zero!
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

if ~exist('method','var') || isempty(method), method='lin'; end
if ~exist('isGraph','var') || isempty(isGraph), isGraph=0; end
if ~exist('sigma','var') || isempty(sigma)
    if strcmpi(method,'ncut')
        sigma=max(D(:))*0.05;
    else
        sigma=max(D(:))*0.5;
    end
end
if ~exist('order','var') || isempty(order), order = 2; end

sparseMode = 0;
if issparse(D)
   sparseMode = 1; 
end

%add this value to nonzero elements
eps=1e-10;
% In case that algorithm changes non-zero element to zero, recover it.
zeroMask = D==0;

switch method
    case 'lin'
        S = 1 - (D./max(D(:))) + eps;
        % Ref:
        % Nascimento, M. & Carvalho, A.
        % A graph clustering algorithm based on a clustering coefficient for weighted graphs
        % Journal of the Brazilian Computer Society, Springer London, 2011, 17, 19-29
    case 'lin2'
        % same as lin, but here the smallest similarity
        % equals the smallest distance after processing.
        minD = min(D(~zeroMask));
        maxD = max(D(:));
        if maxD > 1
            S = D./maxD;
            S = 1 - S + minD/maxD;
        else
            S = 1-D + minD;
        end
        
    case 'ncut'  
        S = exp(-(D./sigma).^order)+eps;
        
    case 'gauss'
        S= exp( - (D.^order ./ (2*sigma^order)) ) +eps;
        
        % reference for gauss2 and gauss3:
        % https://www.aaai.org/Papers/Workshops/2000/WS-00-01/WS00-01-011.pdf
    case 'gauss2'
        S = exp( - (D ./ max( D(:) ) ).^order) +eps;
        
    case 'gauss3'
        S= exp( - D.^order) + eps;
        
    otherwise
        error('dist2sim: wrong mode string!');
end

% set elements to 0 according to zero elements in D (0 indicates no edge
% between two points)
if isGraph
    S(zeroMask) = 0;
end

if sparseMode
    S = sparse(S);
end
