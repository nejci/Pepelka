function S=dist2sim(D,method,sigma,sparseMode)
% For graphs! Zero entries must remain 0 after operation!!!
%
% Converts distance to similarity matrix, using transformation in mode:
% 'lin' - linear, normalized
% 'lin2' - linear, normalized, max=1, min sim. = min. dist
% 'gauss' - gaussian; sigma is half of range if not otherwise defined!
% 'gauss2' - quasi gaussian, normalized, without sigma
% 'gauss3' - quasi gaussian, without sigma

if ~exist('method','var') || isempty(method), method='lin'; end
if ~exist('sigma','var') || isempty(sigma), sigma=max(D(:))*0.5; end
if ~exist('sparseMode','var') || isempty(sparseMode), sparseMode=0; end


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
        
	case 'gauss'		
		S= exp( - (D.^2 ./ (2*sigma^2)) ) +eps;
		
	% reference for gauss2 and gauss3:
	% https://www.aaai.org/Papers/Workshops/2000/WS-00-01/WS00-01-011.pdf
	case 'gauss2'		
		S= exp( - (D ./ max( D(:) ) ).^2) +eps;
		
	case 'gauss3'		
		S= exp( - D.^2) + eps;
		
	otherwise
		error('dist2sim: wrong mode string!');
end

% set elements to 0 according to zero elements in D (0 indicates no edge
% between two points)
S(zeroMask)=0;

if sparseMode
	S=sparse(S);
end
