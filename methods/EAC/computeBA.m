function BA = computeBA(labelsEns)
% BA = computeBA(labelsEns)
% computes binary representation matrix BA of partition ensemble labelsEns
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns   (matrix)	clustering ensemble; each COLUMN corresponds
%							to one clustering (standard in Pepelka).
%
% OUTPUTS
%   BA          (matrix)    binary matrix of size [nClusters X nData], where
%                           BA(i,j)=1 if data point j is a member of cluster i,
%                           and 0 otherwise.
% 
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 11-July-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

PI = labelsEns';

% M = number of clusterings
% N = number of data samples
[M,N] = size(PI);

% save the positions of NaNs
nanmask = isnan(PI);

nClustL = zeros(M,1);
for iM = 1:M
    lbl = PI(iM,:);
    [u,~,iB] = unique(lbl);
    nClustL(iM) = length(u) - sum(nanmask(iM,:));
    
    % repair labels that are not in proper sequential form
    if ~isequal(u, 1:nClustL(iM))
        fprintf(1,'Labels in line %d have been repaired!\n',iM);
        PI(iM,:) = iB;
    end
end
PI(nanmask) = nan;

nClust = sum(nClustL);

% Construct binary representation of PI and save it to P.
% Each cluster is a row with ones on indeces of containing data points.
BA = false(nClust,N);
sep=[0;cumsum(nClustL)];

for p=1:M
    Ptmp = false(nClustL(p),N);
    
    step = nClustL(p);
    linInd = PI(p,:) + (0:step:step*N-1);
    linInd(isnan(linInd))=[];
    Ptmp(linInd) = 1;
    BA( (sep(p)+1) : sep(p+1) ,:) = Ptmp;
end