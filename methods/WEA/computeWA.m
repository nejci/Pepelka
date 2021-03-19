function WA = computeWA(labelsEns,data,dataMode,dataDist,outMode,normalize,weights)
% CO = computeCO(P,numPartitions,inMode,outMode,weights) 
% computes co-occurence accumulation matrix from binary representations or 
% integer labels of clusterings in the ensemble
%--------------------------------------------------------------------------
% INPUTS
% labelsEns     (matrix)  clustering ensemble; each COLUMN corresponds
%                         to the labels of one clustering. Labels must be
%                         positive integers from 1 to k, where k is the
%                         number of clusters.
% data          (matrix)  original data from which labelsEns were generated
% dataMode      (string)  data is:
%                         'data' - data matrix [N X dim]
%                         'dist' - dissimilarity matrix [N X N]
%                         'sim'  - similarity matrix [N X N]
%
% dataDist      (scalar)  if dataMode is 'data', dataDist is used to
%                         calculate pair-wise distances between data. Can
%                         be any distance supported by pdist function. 
%                         If dataMode is 'dist', dataDist is used to
%                         compute similarities.
%
% outMode       - 'full'  : output is N X N matrix (default)
%               - 'vec'   : output is N*(N-1)/2 long vector of pairs%                           
%               - 'sparse': output is sparse matrix
%
% normalize     Should cluster properties be normalized? Default is 1.                            
% weights       optional
%               - vector of weights for each partition if inMode is labels;
%               - if empty, we will use a vector of ones.
%
% OUTPUTS
% WA            weighted co-occurence accumulation values in a form
%               determined by outMode parameter
%
%------- REFERENCE --------------------------------------------------------
% Vega-Pons, S., Ruiz-Shulcloper, J., & Guerra-Gand�n, A. (2011). Weighted
% association based methods for the combination of heterogeneous
% partitions. Pattern Recognition Letters, 32(16), 2163�2170.
% doi:10.1016/j.patrec.2011.05.006
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 13-September-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

[N,ensSize] = size(labelsEns);

if ~exist('normalize','var') || isempty(normalize)
    normalize = 1;
end

% data handling
if ~exist('dataMode','var') || isempty(dataMode)
    dataMode = 'data';
end
if ~exist('dataDist','var') || isempty(dataDist)
    dataDist = 'euclidean';
end

if strcmpi(dataMode,'data')
    data = squareform(pdist(data,dataDist));
    dataMode = 'dist';
end

if strcmpi(dataMode,'dist')
    switch dataDist
        case {'euclidean', 'cityblock', 'minkowski', 'chebychev'}
            data = 1./(1+data);
        case {'correlation', 'spearman'}
            data = 2-data; %correlation coeficient on interval [0,2] instead of [-1,1]
        case {'cosine', 'hamming','jaccard'}
            data = 1-data; %cosine coeficient on interval [0,1]
        otherwise
            error(['Do not know how to calculate similarities from distance ', dataDist])
    end
end

% make diagonal elements 0
data(1:N+1:N^2)=0;

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



% compute number of clusters in each clustering and number of elements in
% each cluster
clustNum = max(labelsEns,[],2);
clustSizes = zeros(ensSize,max(clustNum));
for c=1:ensSize
    clustSizes(c,1:clustNum(c)) = histc(labelsEns(c,:),1:clustNum(c));
end

if normalize
    clustNum = clustNum ./ max(clustNum);
    data = data./max(data(:));
    clustSizes = clustSizes ./ max(clustSizes(:));
end

switch outMode
    case 'full'
        WA = zeros(N,N);
                
        for i=1:N-1
            m = bsxfun(@eq,labelsEns,labelsEns(:,i));
            cSize = clustSizes(sub2ind(size(clustSizes),(1:ensSize)',labelsEns(:,i)));            
            d = data(i,:);
            m = bsxfun(@times,m,clustNum);
            m = bsxfun(@times,m,d);
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            WA(i,:)= sum(m,1);
        end
        WA(end,:) = WA(:,end)';
        %WA = triu(WA,1); % upper part of the matrix without diagonal
        %WA = WA + WA' + eye(N,N);
        
    case 'vec'
        nPairs = N*(N-1)/2;
        WA = zeros(1,nPairs);
        ind=[0 cumsum(N-1:-1:1)];
        
        for i=1:N-1
            m = bsxfun(@eq,labelsEns,labelsEns(:,i));
            cSize = clustSizes(sub2ind(size(clustSizes),(1:ensSize)',labelsEns(:,i)));            
            d = data(i,:);
            m = bsxfun(@times,m,clustNum);
            m = bsxfun(@times,m,d);
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            s = sum(m,1);
            WA( ind(i)+1 : ind(i+1)) = s(i+1:end);
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
            d = data(i,:);
            m = bsxfun(@times,m,clustNum);
            m = bsxfun(@times,m,d);
            m = bsxfun(@rdivide,m,cSize);
            m = bsxfun(@times,m,weights);
            s = sum(m,1);
            AMvec( ind(i)+1 : ind(i+1)) = s(i+1:end);
        end
        WA = sparse(I,J,AMvec,N,N);
        WA = WA + WA';
        
    otherwise
        error('Wrong outMode value!');
end