function GDI = indexGDI(data,labels,interDist,intraDist,distM,centers,D2C,C2C)

% GDI = INDEXGDI(data,labels,graph_type,options)
%--------------------------------------------------------------------------
% Generalized Dunn cluster validity index.
% Best cluster partition maximizes the index value.
%--------------------------------------------------------------------------
% INPUTS
%   data        (matrix)    data matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%
%	interDist   (string/int) distance between two clusters
%                               1 - 'single' (original Dunn)
%                               2 - 'complete'
%                               3 - 'average'
%                               4 - 'centroid'
%                               5 - 'avg2cent'
%                               6 - 'hausdorff'
%                               - if empty, compute all of them
%
%   intraDist   (string/int) cluster diameter:
%                               1 - 'complete' (original Dunn)
%                               2 - 'average'
%                               3 - 'avg2cent'
%                               - if empty, compute all of them
%
%   distM  		(matrix)	dissimilarity matrix [n X n] OR
%               (string)    an distance identifier that is passed to pdist
%                           function (default: 'euclidean').
%   centers     (matrix)    centers of clusters [K X dim]; if empty, means
%                           of data points in cluster are its center
%   D2C         (vector)    distance from each data point to its cluster
%                           center; if empty, Euclidean distance is
%                           computed.
%               (string)    an distance identifier that is passed to pdist2
%                           function (default: 'euclidean').
%   C2C         (matrix)    pairwise distances between cluster centers; 
%                           if empty, Euclidean distance is computed.
%               (string)    an distance identifier that is passed to pdist
%                           function (default: 'euclidean').
%--------------------------------------------------------------------------
% OUTPUTS:
%   GDI         (scalar)	value of index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Bezdek, J. C., & Pal, N. R. (1998). Some new indexes of cluster validity.
% IEEE Transactions on Systems, Man, and Cybernetics, 28 , 301�315.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 28-May-2014 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

[N,dim] = size(data);
dtype = 'euclidean';
distDefined = 0;

% compute dissimilarity matrix, if not given as input argument
if ~exist('distM','var') || isempty(distM)
    distM = squareform(pdist(data,'euclidean'));
    distDefined=1;
end
if ischar(distM)
    dtype = distM;
    distDefined=1;
    distM = squareform(pdist(data,distM));   
end

K = max(labels);

% determine which inter-cluster distances to compute
interAll = {'single','complete','average','centroid','avg2cent','hausdorff'};
interAllNum = 1:6;
if ~exist('interDist','var') || isempty(interDist)
    interSel = [1 1 1 1 1 1];
    interPos = 1:6;
    interNum = 6;
else
    if isnumeric(interDist)
        [interSel,interPos] = ismember(interAllNum, interDist);
    else
        if ~iscellstr(interDist)
            interDist = {interDist};
        end
        [interSel,interPos] = ismember(interAll, interDist);
    end
    interNum = length(interDist);
end


% determine which intra-cluster distances (diameters) to compute
intraAll = {'complete','average','avg2cent'};
intraAllNum = 1:3;
if ~exist('intraDist','var') || isempty(intraDist)
    intraSel = [1 1 1];
    intraPos = 1:3;
    intraNum = 3;
else
    if isnumeric(intraDist)
        [intraSel,intraPos] = ismember(intraAllNum, intraDist);
    else
        if ~iscellstr(intraDist)
            intraDist = {intraDist};
        end
        [intraSel,intraPos] = ismember(intraAll, intraDist);
    end
    intraNum = length(intraDist);
end



% Do we need cluster centers, C2C and D2C?
if interSel(4) || interSel(5) || intraSel(3)
    if ~exist('centers','var') || isempty(centers)
        % compute centers of clusters
        centers = zeros(K,dim);
        for i = 1:K
            Dtmp = data(labels==i,:);
            centers(i,:) = mean(Dtmp,1);
        end
    end
    % compute distance from data to each cluster center
    if (interSel(5) || intraSel(3))  
        if (~exist('D2C','var') || isempty(D2C))
            D2C = dtype;
        end
        if ischar(D2C)
            if distDefined && ~strcmpi(D2C,dtype) 
                error('Different distances for data and clusters!');
            end
            D2C_dist = D2C;
            D2C = zeros(N,K);
            
            for i = 1:K            
                D2C(:,i) = pdist2(data,centers(i,:),D2C_dist);
            end
        end
    end
    % compute pairwise distances between cluster centers
    if interSel(4)
        if ~exist('C2C','var') || isempty(C2C)
            C2C = dtype;
        end
        if ischar(C2C)
            if distDefined && ~strcmpi(C2C,dtype)
                error('Different distances for data and clusters!');
            end
            C2C = squareform(pdist(centers,dtype));
        end
    end
end


interMat = zeros(interNum,(K^2-K)/2);
intraMat = zeros(intraNum,K);

kInd = 1;
for k1 = 1:K
    lab1 = labels == k1;
    
    % compute cluster diameter
    Dtmp = distM(lab1,lab1);
    nk1 = size(Dtmp,1);
    
    if intraSel(1)
        intraMat(intraPos(1),k1) = max(Dtmp(:));
    end
    if intraSel(2)
        if nk1 == 1
            intraMat(intraPos(2),k1) = 0;
        else
            intraMat(intraPos(2),k1) = sum(Dtmp(:)) / (nk1^2-nk1);
        end
    end
    if intraSel(3)
        intraMat(intraPos(3),k1) = 2*(sum(D2C(lab1,k1)) / nk1);
    end
    
    % compute distances between pairs of clusters
    for k2 = (k1+1):K
        lab2 = labels == k2;
        Dtmp = distM(lab1,lab2);
        nk2 = size(Dtmp,2);
        
        if interSel(1)
            interMat(interPos(1),kInd) = min(Dtmp(:));
        end
        if interSel(2)
            interMat(interPos(2),kInd) = max(Dtmp(:));
        end        
        if interSel(3)
            interMat(interPos(3),kInd) = sum(Dtmp(:)) / (nk1*nk2);
        end
        if interSel(4)
            interMat(interPos(4),kInd) = C2C(k1,k2);
        end
        if interSel(5)
            interMat(interPos(5),kInd) = (sum(D2C(lab1,k2)) + sum(D2C(lab2,k1))) / (nk1+nk2);
        end
        if interSel(6)
            % Compute Hausdorff distance between two sets of points  
            % Based on the code of Zachary Danziger (FEX 26738)
            % Obtain the value of the point, p, in P with the largest minimum distance
            % to any point in Q.
            vp = max(min(Dtmp,[],2));
            % Obtain the value of the point, q, in Q with the largest minimum distance
            % to any point in P.
            vq = max(min(Dtmp,[],1));
            % choose maximum value between vp and vq
            hausdorff = max(vp,vq);            
            interMat(interPos(6),kInd) = hausdorff;            
        end
        
        kInd = kInd+1;
    end
end

GDI = bsxfun(@rdivide, min(interMat,[],2), max(intraMat,[],2)');


% denominator=[];
% 
% for i2=1:K
%     indi=find(labels==i2);
%     indj=find(labels~=i2);
%     x=indi;
%     y=indj;
%     temp=distM(x,y);
%     denominator=[denominator;temp(:)];
% end
% 
% num=min(min(denominator));
% neg_obs=zeros(size(distM,1),size(distM,2));
% 
% for ix=1:K
%     indxs=find(labels==ix);
%     neg_obs(indxs,indxs)=1;
% end
% 
% dem=neg_obs.*distM;
% dem=max(max(dem));
% 
% DI=num/dem;
end