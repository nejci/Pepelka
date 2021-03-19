function [APN, AD, ADM, FOM] = stability(data, distMat, labels, labelsDel, distMethod)

% [APN, AD, ADM, FOM] = stability(data, distMat, labels, labelsDel, distMethod)
%--------------------------------------------------------------------------
% Copied from description of clValid:
% "The stability measures are a special version of internal measures which 
% evaluate the stability of a clustering result by comparing it with the 
% clusters obtained by removing one column at a time. These measures include 
% the average proportion of non-overlap (APN), the average distance (AD), 
% the average distance between means (ADM), and the figure of merit (FOM). 
% The APN, AD, and ADM are all based on the cross-classification table of 
% the original clustering with the clustering based on the removal of one 
% column. The APN measures the average proportion of observations not placed 
% in the same cluster under both cases, while the AD measures the average 
% distance between observations placed in the same cluster under both cases 
% and the ADM measures the average distance between cluster centers for 
% observations placed in the same cluster under both cases. 
% The FOM measures the average intra-cluster variance of the deleted column, 
% where the clustering is based on the remaining (undeleted) columns. 
% 
% In all cases the average is taken over all the deleted columns, and all 
% measures should be minimized."
%--------------------------------------------------------------------------
% INPUTS
%   data			(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels			(vector)	array of non-negative integers determining
%								the labels of data samples when clustering
%								whole data
%
%	labelsDel		(matrix)	matrix [n X d] of labels of data when 
%								clustering data without corresponding
%								column, i.e. labelsDel(:,i) contains labels
%								when clustering data without column i.
%
%	distMethod		(string)	which type of distance to use. See
%                               help('pdist') for the list of supported 
%                               metrics.
%								
%--------------------------------------------------------------------------
% OUTPUTS:
%   value			(vector)	value of indeces APN, AD, ADM, and FOM
%
%--------------------------------------------------------------------------
% REQUIRES:         Statistics Toolbox (crosstab, nanmean)
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
% Part of Pepelka package based on R package clValid.
%
%------- REFERENCE --------------------------------------------------------
% Brock, G., Pihur, V., Datta, S. and Datta, S. (2008), 
% clValid: An R Package for Cluster Validation, 
% Journal of Statistical Software 25(4) http://www.jstatsoft.org/v25/i04
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 4-Oct-2012 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================

if (~all(size(data) == size(labelsDel)))
    error('data and labelsDel must be of the same size!');
end

% compute distance matrix of data if none provided
if (isempty(distMat))
    distMat = squareform(pdist(data,distMethod));
end

anyNan = any(any(isnan(data)));

[n,d] = size(data);

stabmeas = zeros(4,d);

% loop over deleted columns
c1 = 1:max(labels);
nc1 = length(c1);

for del = 1:d
    labelsDel_i = labelsDel(:,del);
    
    c2 = unique(labelsDel_i)';    
    nc2 = length(c2);
    
    overlap = crosstab(labels,labelsDel_i);
    
    dij = zeros(nc1,nc2);
    dij2 = zeros(nc1,nc2);
    
    ii = 1;
    
    for i = c1
        jj = 1;
        
        ci = (labels == i);
        xbari = nanmean(data(ci,:),1);
        
        for j = c2
           
           cDelj = (labelsDel_i == j);
           
           cl = sum(ci)*sum(cDelj);
           
           if (cl > 0)
               dij(ii,jj) = mean(nanmean(distMat(ci, cDelj) ));
           end
           xbarj = nanmean(data(cDelj,:),1);
           differ = xbari - xbarj;
           
           if (~isempty(differ))
               if(anyNan)
                   differ = differ(~isnan(differ));
                   dij2(ii,jj) = sqrt(mean(differ.^2));
               else
                   dij2(ii,jj) = sqrt(sum(differ.^2));
               end
           else
               dij2(ii,jj) = 0;
           end
           jj = jj + 1;
        end
        
        ii = ii + 1;                            
    end
    
    rs = repmat(sum(overlap,2),1,nc2);
    cs = repmat(sum(overlap,1),nc1,1);
    
    xbar_vec = zeros(n,1);
    for k = c2
       xbar = mean(data(labelsDel_i == k, del));
       xbar_vec(labelsDel_i == k) = xbar;
    end
    
    APNi = 1 - sum(sum(overlap.^2 ./ rs)) /sum(overlap(:));
    ADi = sum(overlap(:) .* dij(:)) / n;
    ADMi = sum(overlap(:) .* dij2(:)) / n;
    FOMi = sqrt(nanmean((data(:,del) - xbar_vec).^2))/sqrt((n - nc1)/n);
    
    stabmeas(:,del) = [APNi,ADi,ADMi,FOMi];
end

val = mean(stabmeas,2);
APN = val(1);
AD = val(2);
ADM = val(3);
FOM = val(4);
