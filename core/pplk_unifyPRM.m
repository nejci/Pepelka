function PRM = pplk_unifyPRM(V, mask, method)
% PRM = pplk_unifyPRM(V, mask, method)
%
% V - matrix of indices' values
% mask - logic vector with 1's where reverse of index's value has to be
% considered
% method -	
%   'minmax': max-like: value_i/max, min-like: min/value_i
%   'range' : normalize on [0,1] (default)
%   'range0': move range to start with 0
%   'range1': move range to end with 1
%   'pos'   : move indeces with negative values to start with 0
%   'prob'  : transform to probabilities, sum is 1
%   'rank'	: rank values, 1 for the poorest value, 2, 3, etc.
%   'rank10': rank values on interval [0,1]; 1 for the highest value
%   'rankRankAggreg': rank values, 1 for the poorest value, other are less
%             than 1 (suitable for Rank Aggregation)
%   'none'  : do nothing, PRM = V
%


% If not specifed otherwise, we assume that every index assigns low value
% to low quality clustering and high to better ones.
if ~exist('mask','var') || isempty(mask)
	mask = false(1,size(V,2));
end

if length(mask) ~= size(V,2)
   error('Length of mask vector must equal the number of indeces.'); 
end

if ~exist('method','var') || isempty(method)
	method = 'range';
end

if strcmpi(method,'none')
   PRM = V;
   return;
end


ensembleSize = size(V,1);

switch method
    case 'minmax'        
        %   if min value of a column is negative, shift all values by
        %   this min value
        mask_neg = min(V,[],1) < 0;
        tmp = V(:,mask_neg);
        tmp = bsxfun(@plus, tmp, -min(tmp,[],1));
        V(:,mask_neg) = tmp;
        
        % max-like: value_i/max
        PRM = bsxfun(@rdivide,V,max(V,[],1));
        
        % min-like: min/value_i: be careful with 0 (0/x = 0)
        % shift all min-like columns with 0's by 1
        mask_zero = (min(V,[],1) == 0) & mask;
        tmp = V(:,mask_zero);
        tmp = bsxfun(@plus, tmp, 1);
        V(:,mask_zero) = tmp;
               
        PRM(:,mask) = bsxfun(@rdivide,min(V(:,mask),[],1),V(:,mask));
        
        % special case, when a column has a constant value
        mask_constant = all(bsxfun(@eq,V,V(1,:)),1);
        PRM(:,mask_constant) = 1;
        
	case 'rangeOld' % OLD: does not scale to [0,1]
        %   if min value of a column is negative, shift all values by
        %   this min value
        mask_neg = min(V,[],1) < 0;
        tmp = V(:,mask_neg);
        tmp = bsxfun(@plus, tmp, -min(tmp,[],1));
        V(:,mask_neg) = tmp;
        
		% normalize each column on [0,1]:   (x - min(x))./ max(x)
		PRM = bsxfun(@minus,V,min(V,[],1));
		PRM = bsxfun(@rdivide,PRM,max(V,[],1));
		PRM(:,mask) = 1-PRM(:,mask);
        
        % special case, when a column has a constant value
        mask_constant = all(bsxfun(@eq,V,V(1,:)),1);
        PRM(:,mask_constant) = 1;
	
    case 'range' 
        %   if min value of a column is negative, shift all values by
        %   this min value
        mask_neg = min(V,[],1) < 0;
        tmp = V(:,mask_neg);
        tmp = bsxfun(@plus, tmp, -min(tmp,[],1));
        V(:,mask_neg) = tmp;
        
		% normalize each column on [0,1]:   (x - min(x))./ (max(x) -
		% min(x))
		PRM = bsxfun(@minus,V,min(V,[],1));
		PRM = bsxfun(@rdivide,PRM,max(PRM,[],1));
		PRM(:,mask) = 1-PRM(:,mask);
        
        % special case, when a column has a constant value
        mask_constant = all(bsxfun(@eq,V,V(1,:)),1);
        PRM(:,mask_constant) = 1;
        
	case 'range0'
		% move each column on [0,maxV-minV]
		PRM = bsxfun(@minus,V,min(V,[],1));
		PRM(:,mask) = bsxfun(@minus, max(PRM(:,mask),[],1), PRM(:,mask));
		
	case 'range1'
		% move each column on [minV-maxV+1,1]
		d = bsxfun(@minus,max(V,[],1),1);
        PRM = bsxfun(@minus, V, d);
		PRM(:,mask) = bsxfun(@plus, 1-PRM(:,mask) ,min(PRM(:,mask),[],1));
		
    case 'pos'
        %   bound lower values to 0; if min value of column is negative, shift
        %   all values by this min value
        mask_neg = min(V,[],1) < 0;
        tmp = V(:,mask_neg);
        tmp = bsxfun(@plus, tmp, -min(tmp,[],1));
        V(:,mask_neg) = tmp;
        
        % mirror masked indices over mid-point
        tmp=V(:,mask);
        mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
        tmp = tmp + 2 * bsxfun(@minus, mid, tmp);
        V(:,mask) = tmp;
        
        PRM = V;
		
	case 'prob'
        %   bound lower values to 0; if min value of column is negative, shift
        %   all values by this min value
        mask_neg = min(V,[],1) < 0;
        tmp = V(:,mask_neg);
        tmp = bsxfun(@plus, tmp, -min(tmp,[],1));
        V(:,mask_neg) = tmp;        
        
        % mirror masked indices over mid-point
        tmp=V(:,mask);
        mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
        tmp = tmp + 2 * bsxfun(@minus, mid, tmp);
        V(:,mask) = tmp;
        
        % calculate probabilites: sum of columns has to be 1
		PRM = bsxfun(@rdivide, V, sum(V,1));
        
        % special case, when a column has a constant value
        mask_constant = all(bsxfun(@eq,V,V(1,:)),1);
        PRM(:,mask_constant) = 1/ensembleSize;
		
		
	case 'rank'
		% rank each column value (index score) with number from 1 to
		% ensembleSize (number of observations) according to its score.
        % Better score = higher rank ("inverse rank").
		% Mind repeated values - calculate average rank.
		% Consider the behaviour of certain index (min/max value means
		% what).
		% Example: Selected indices: DB (min-like), DN (max-like). 
		% PRM = [32, 55, 0.2, 12; 0.3 0.7 5.6 5.6]'
		% ranks (numbers 1:4): [2, 1, 4, 3; 1, 2, 3.5, 3.5];
		
		% function tiedrank ranks the smallest number with 1. The higher
		% the index score, the higher the rank gets. It is ok for MAX-like
		% indices. For MIN-like ones, we have to reverse ordering 
		% (subtract rank from ensembleSize+1).
		PRM = tiedrank(V);
		PRM(:,mask) = (ensembleSize+1) - PRM(:,mask);
        
	case 'rank10'
		% rank each column value (index score) with number from 1/ensembleSize (worst)
		% to 1 (best) according to its score.
        % Better score = higher rank ("inverse rank").
		% Mind repeated values - calculate average rank.
				
		% function tiedrank ranks the smallest number with 1. The higher
		% the index score, the higher the rank gets. It is ok for MAX-like
		% indices. For MIN-like ones, we have to reverse ordering 
		% (subtract rank from ensembleSize+1).
		PRM = tiedrank(V);
		PRM(:,mask) = (ensembleSize+1) - PRM(:,mask);
        PRM = PRM /ensembleSize; % scale on [0,1]
        
    case 'rankRankAggreg'
		% Convenient for robustRankAggregation
        % same as rank10 but with inverse notion:
        % rank each column value (index score) with number from 0 (best)
		% to 1 (worst) according to its score.
        % NOTE: smaller rank means better result!
		PRM = tiedrank(V);
        PRM(:,mask) = PRM(:,mask) / ensembleSize;
		PRM(:,~mask) = ((ensembleSize+1) - PRM(:,~mask)) / ensembleSize;
        
	otherwise
		error('Wrong unification method!');
end
end
