function val = indexI(data, labels, p, V)
% val = indexI(data, labels, p)
%--------------------------------------------------------------------------
% Index I: (internal) cluster validity index
% Best parition maximizes the index value.
%--------------------------------------------------------------------------
% INPUTS
%   data			(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels			(vector)	array of non-negative integers determining
%								the labels of data samples
%					(matrix)	fuzzy member matrix [k X n]
%
%	p				(scalar)	power used to control the contrast between
%								the different cluster configurations 
%                               (default p=2).
%   V               (matrix)    cluster centers [k X d]
%--------------------------------------------------------------------------
% OUTPUTS:
%   val				(scalar)	value of index I
%
%------- CALCULATION ------------------------------------------------------
% Index I is defined as
% I = (1/k * E1/Ek * Dk)^p  ,
% where
% k is number of clusters
% E1 = sum_{i=1}^n(||data_i - v_D||)  (sum of Euclidean distances between each data sample and data mean)
% Ek = sum_{i=1}^k(sum_{j = 1}^n ( U(i,j) ||x_j - v_i|| ))  (sum of Euclidean distances (weighted by membership matrix) between each data sample and each cluster center)
% Dk = max_{v_i,v_j \in V} x(||v_i - v_j||)   (maximum distance between cluster centers)
% p is parameter, default is 2
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Maulik, U., & Bandyopadhyay, S. (2002).
% Performance Evaluation of Some Clustering Algorithms and Validity Indices.
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% 24 , 1650-1654.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 22-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[n,d]=size(data);
k = max(labels);

if isvector(labels)
	E=eye(k);
	U = E(:,labels);
else
	U = labels; % Fuzzy membership matrix U
end

if ~exist('p','var') || isempty(p)
	p = 2;
end

if ~exist('V','var') || isempty(V)
    % find cluster means or centers
    V = zeros(k,d);
    for c=1:k
        V(c,:)=mean(data(labels==c,:),1);
    end
end

% sum of distances between data and cluster centers
Ek=0;
for i=1:k    
    dsum = U(i,:) * sqrt(sum(bsxfun(@minus,data,V(i,:)).^2,2));
    Ek = Ek + dsum;
end

disttc = sqrt(dist_euclidean(V,V));
Dk = max(max(disttc));

% find center (mean) of the entire data
vD = mean(data,1);

% compute sum of distances between data mean and data points 
E1 = sum(sqrt(sum(bsxfun(@minus,data,vD).^2,2)));

% final score
val =((E1*Dk)/(k*Ek))^p;
