function [XB,XBmod] = indexXB(data,labels)

% val = indexXB(data, labels, mod)
%--------------------------------------------------------------------------
% Index XB: (internal) cluster validity index
% Best parition minimizes the index value.
%--------------------------------------------------------------------------
% INPUTS
%   data			(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels			(vector)	array of non-negative integers determining
%								the labels of data samples
%					(struct)	.U: fuzzy membership matrix [k X n]
%                               .center: cluster centers [k X d]
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   XB				(scalar)	value of index XB
%   XBmod			(scalar)	value of modified index XB
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
% Part of Pepelka package.
%
%------- REFERENCES -------------------------------------------------------
% Xie, X. L., & Beni, G. (1991). A validity measure for fuzzy clustering.
% Pattern Analysis and Machine Intelligence, IEEE Transactions on, 13(8),841-847
%
% Kim, M., & Ramakrishna, R. S. (2005). New indices for cluster validity
% assessment. Pattern Recognition Letters, 26(15), 2353-2363.
%------- VERSION ----------------------------------------------------------
% Version: 1.2
% Last modified: 25-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================


[n,d]=size(data);


if isstruct(labels)
    clusterCenters = labels.center;
    U = labels.U;
    k=size(U,1);
    
else
    if isvector(labels)
        k = max(labels);
        E=eye(k);
        U = E(:,labels);
        clsize = sum(U,2);      
        clusterCenters = bsxfun(@rdivide, U*data, clsize);  
        % find cluster means/centers
%         clusterCenters = zeros(k,d);
%         for c=1:k
%             ind = logical(U(c,:));
%             clusterCenters(c,:)=mean(data(ind,:),1);
%         end
        
    else
        error('Variable labels has to be struct or vector.');
    end
end


% distances (squared Euclidean) between data samples and centers
distCenter2Data = distSqEucl(clusterCenters,data);

j_XB = sum(sum(distCenter2Data .* (U.^2),2))/n;

ni = sum(U>0, 2);
j_XBmod = max( sum(distCenter2Data .* (U.^2),2) ./ ni );


% Squared Euclidean distances between cluster centers
distCenter2Center=distSqEucl(clusterCenters,clusterCenters);
distCenter2Center(1:k+1:end)=Inf; % diagonal elements are not interesting

% minimal distance between cluster centers
sep = min(min(distCenter2Center));

XB    = j_XB/sep;
XBmod = j_XBmod/sep;

end

function out = distSqEucl(center, data)
%	OUT = distSqEucl(CENTER, DATA) calculates the Squared Euclidean distance
%	between each row in CENTER and each row in DATA, and returns a
%	distance matrix OUT of size M by N, where M and N are row
%	dimensions of CENTER and DATA, respectively, and OUT(I, J) is
%	the distance between CENTER(I,:) and DATA(J,:).

%	Roger Jang, 11-22-94, 6-27-95.
%       Copyright 1994-2002 The MathWorks, Inc.
%       $Revision: 1.13 $  $Date: 2002/04/14 22:20:29 $
%  Modified by Nejc Ilc

out = zeros(size(center, 1), size(data, 1));

% fill the output matrix

if size(center, 2) > 1,
    for k = 1:size(center, 1),
        out(k, :) = sum(((data-ones(size(data, 1), 1)*center(k, :)).^2),2)';
    end
else	% 1-D data
    for k = 1:size(center, 1),
        out(k, :) = abs(center(k)-data)';
    end
end
end
