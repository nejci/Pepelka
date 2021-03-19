function [newE, no_allcl] = relabelCl(E) 
%==========================================================================
% FUNCTION: [newE, no_allcl] = relabelCl(E)
% DESCRIPTION:	This function is used for relabelling clusters in the
% ensemble 'E'.
% %stevilke gruc v 2. in vseh nadaljnih razvrstitvah povecamo za stevilo
% gruc v prejsnji razvrstitvi. 
%				
%
% INPUTS:    E = N-by-M matrix of cluster ensemble
%
% OUTPUT: newE = N-by-M matrix of relabeled ensemble
%     no_allcl = total number of clusters in the ensemble
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

[N,M] = size(E); % no. of clustering
newE = zeros(N,M);

%--- first clustering
ucl = unique(E(:,1)); % all clusters in i-th clustering
if (max(E(:,1) ~= length(ucl))) %NAPAKA??? Prav: max(E(:,1)) ~= length(ucl) ?
    for j = 1:length(ucl)
        newE(E(:,1) == ucl(j),1) = j; % re-labelling
    end
end

%--- the rest of the clustering
for i = 2:M
    ucl = unique(E(:,i)); % all clusters in i-th clustering
    prevCl = length(unique(newE(:,[1:i-1])));
    for j = 1:length(ucl)
        newE(E(:,i) == ucl(j),i) = prevCl + j; % re-labelling (add an offset)
    end
end

no_allcl = max(max(newE));