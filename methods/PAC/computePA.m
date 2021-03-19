function PA = computePA(labelsEns,dim,outMode,weights)
% PA = computePA(labelsEns,dim,outMode,weights)
% computes probability co-occurence accumulation matrix from  
% integer labels of clusterings in the ensemble
%--------------------------------------------------------------------------
% INPUTS
% labelsEns     (matrix)  clustering ensemble; each COLUMN corresponds
%                         to the labels of one clustering. Labels must be
%                         positive integers from 1 to k, where k is the
%                         number of clusters.
% dim           (scalar)  dimensionality of the clustered data
%
% outMode       - 'full'  : output is N X N matrix (default)
%               - 'vec'   : output is N*(N-1)/2 long vector of pairs
%                           (diagonal elements (ones) are not included)
%               - 'sparse': output is sparse matrix
%                           (diagonal elements (ones) are not included)
% weights       optional
%               - vector of weights for each partition if inMode is labels;
%               - if empty, we will use a vector of ones.
%
% OUTPUTS
% PA            co-occurence probability accumulation values in a form
%               determined by outMode parameter
%
%------- REFERENCE --------------------------------------------------------
% Wang, X., Yang, C., & Zhou, J. (2009). Clustering aggregation by
% probability accumulation. Pattern Recognition, 42, 668�675.
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 9-September-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

[N,ensSize] = size(labelsEns);

if ~exist('outMode','var') || isempty(outMode)
    outMode = 'full';
end

if ~exist('weights','var') || isempty(weights)
    weights = ones(ensSize,1);
end

if any(isnan(labelsEns(:)))
    error('Labels cannot be NaN!');
end

labelsEns = labelsEns';

% compute size of clusters in each clustering
clustMax = max(labelsEns,[],2);
clustSizes = zeros(ensSize,max(clustMax));
for c=1:ensSize
    clustSizes(c,1:clustMax(c)) = histc(labelsEns(c,:),1:clustMax(c));
end


switch outMode
    case 'full'
        PA = zeros(N,N);
                
        for i=1:N-1
            m = bsxfun(@eq,labelsEns,labelsEns(:,i));
            cSize = clustSizes(sub2ind(size(clustSizes),(1:ensSize)',labelsEns(:,i)));
            cSize = 1+ cSize.^(1/dim);
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            PA(i,:)= sum(m,1);
        end
        PA = triu(PA,1); % upper part of the matrix without diagonal
        PA = PA + PA' + eye(N,N);
        
    case 'vec'
        nPairs = N*(N-1)/2;
        PA = zeros(1,nPairs);
        ind=[0 cumsum(N-1:-1:1)];
        
        for i=1:N-1
            m = bsxfun(@eq,labelsEns,labelsEns(:,i));
            cSize = clustSizes(sub2ind(size(clustSizes),(1:ensSize)',labelsEns(:,i)));
            cSize = cSize.^(1/dim) +1;
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            s = sum(m,1);
            PA( ind(i)+1 : ind(i+1)) = s(i+1:end);
        end        
        
    case 'sparse'
        % try if sparse matrix would fit into memory
        AMvec = zeros(1,N*(N-1)/2);
        V = nchoose2(1:N);
        I = V(:,1)';
        J = V(:,2)';
        ind=[0 cumsum(N-1:-1:1)];
        
        for i=1:N-1
            m = bsxfun(@eq,labelsEns,labelsEns(:,i));
            cSize = clustSizes(sub2ind(size(clustSizes),(1:ensSize)',labelsEns(:,i)));
            cSize = cSize.^(1/dim) +1;
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            s = sum(m,1);
            AMvec( ind(i)+1 : ind(i+1)) = s(i+1:end);
        end
        PA = sparse(I,J,AMvec,N,N);
        PA = PA + PA';
        
    otherwise
        error('Wrong outMode value!');
end