function [SSI, SSIcum]= indexSSI(D,labels,centersMode)
% [SSI, SIcum] = INDEXSSI(D,labels,centersMode)
%--------------------------------------------------------------------------
% Cluster internal validity index Simple Structure Index (Dolnicar et al., 1999).
%--------------------------------------------------------------------------
% INPUTS
%   D       	(matrix)	matrix [N X D] with N D-dimensional samples OR
%                           matrix with K D-dimensional centers
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   centersMode (bool)      0 - data is data matrix (default)
%                           1 - data is centers matrix
%--------------------------------------------------------------------------
% OUTPUTS:
%   SSI         (scalar)	value of the SSI index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Dolnicar, S., Grabler, K., & Mazanec, J. (1999). A tale of three cities:
% Perceptual charting for analysing destination images. In A. Woodside
% (Ed.), Consumer psychology of tourism, hospitality and leisure (pp.
% 39-62). London, U.K.: CAB International.
%
% Dimitriadou, E., Dolnicar, S., & Weingessel, A. (2002). An examination of
% indexes for determining the number of clusters in binary data sets.
% Psychometrika, 67(1), 137�159. doi:10.1007/BF02294713
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 8-8-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

K = max(labels);

% compute binary meembership matrix U
E = eye(K);
U = logical(E(:,labels));
clsize = sum(U,2)';

% compute centers if not specified
if ~exist('centersMode','var') || isempty(centersMode)
    centersMode = 0;
end

if ~centersMode
    % D is data matrix
    centers = bsxfun(@rdivide, U*D, clsize');    
else
    % D is centers matrix
    centers = D;
end


[ncl,nvar] = size(centers);
cmax = max(centers,[],1);
cmin = min(centers,[],1);
[~,cord] = sort(centers,1);
cmaxi = cord(ncl,:);
cmini = cord(1,:);

meanmean = mean(centers(:));
absmdif = abs(mean(centers,1) - meanmean);
span = cmax - cmin;
csizemax = clsize(cmaxi);
csizemin = clsize(cmini);

hiest = nvar;
hiestw = hiest * max(max(csizemax), max(csizemin)) * exp(-min(absmdif));

SSIcum = sum(span)/hiest;
SSI = (span .* exp(-absmdif)) * (sqrt(csizemax.*csizemin) / hiestw)';


