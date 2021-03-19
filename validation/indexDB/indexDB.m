function val = indexDB(data,labels,mod,p,q,V)
% val = indexDB(data,labels,mod,p,q)
%--------------------------------------------------------------------------
% Cluster internal validity index - Davies-Bouldin index
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   mod         (scalar)    if 0 run default algorithm, otherwise
%                           the modification (Kim & Ramakrishna, 2005).
%   p           (scalar)    parameter of the Minkowski distance
%                           p = 1: cityblock distance
%                           p = 2: Euclidean distance (default)
%                           ...
%   q           (scalar)    parameter of dispersion measure of a cluster
%                           q = 1: average Euclidean measure (default)
%                           q = 2: standard deviation of the distance
%                           ...
%   V           (matrix)    cluster centers [k X d]
%--------------------------------------------------------------------------
% OUTPUTS:
%   val          (scalar)	value of the Davies-Bouldin index
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Davies, D. L., & Bouldin, D. W. (1979). A Cluster Separation Measure.
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 1(2),
% 224�227.
%
% Kim, M., & Ramakrishna, R. S. (2005). New indices for cluster validity
% assessment. Pattern Recognition Letters, 26(15), 2353�2363.
% doi:10.1016/j.patrec.2005.04.007
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 22-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('mod','var') || isempty(mod)
    mod=0;
end
if ~exist('p','var') || isempty(p)
    p=2;
end
if ~exist('q','var') || isempty(q)
    q=1;
end

[n,d]=size(data);
k=max(labels);

E = eye(k);
U = logical(E(:,labels));

% number of data points in clusters
nk = sum(U,2);

if ~exist('V','var') || isempty(V)
    % find cluster means or centers
    V = bsxfun(@rdivide, U*data, nk);
end

% distances between data points and their cluster centers
S=zeros(1,k);
for i=1:k
    dataK = data(U(i,:),:);
    S(i) = (sum(sqrt(sum(bsxfun(@minus,dataK,V(i,:)).^2,2)).^q)/nk(i))^(1/q);
end

% inter cluster distances (between centers)
M = squareform(pdist(V,'minkowski',p));

% deafult algorithm [Davies-Bouldin, 1979]
Rm = zeros(k);
if mod == 0
    for i=1:k-1
        for j=i+1:k
            R = (S(i)+S(j))/M(i,j);
            Rm(i,j) = R;
            Rm(j,i) = R;
        end
    end
    sumR = sum(max(Rm,[],2));
    val = sumR/k;

% modified DB index [Kim & Ramakrishna, 2005]
else
    sumR = 0;
    for i=1:k
        maxSij=0;
        minD = Inf;        
        for j=1:k
            if i~=j
                Sij = S(i)+S(j);
                d = M(i,j);                
                if Sij > maxSij
                    maxSij=Sij;
                end                
                if d < minD && d ~=0 % fix: d must not equal 0
                    minD = d;
                end
            end
        end
        sumR = sumR + maxSij/minD;
    end
    val = sumR/k;
end