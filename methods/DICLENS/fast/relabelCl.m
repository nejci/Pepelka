function [newE, no_allcl, C] = relabelCl(E) 
%==========================================================================
% FUNCTION: [newE, no_allcl] = relabelCl(E)
% DESCRIPTION:	This function is used for relabelling clusters in the
% ensemble 'E'.
%				
%
% INPUTS:    E = N-by-M matrix of cluster ensemble
%
% OUTPUT: newE = N-by-M matrix of relabeled ensemble
%     no_allcl = total number of clusters in the ensemble
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

newE = E;
nClsOrig = max(newE,[],1);
C = [0 cumsum(nClsOrig)]; 
newE = bsxfun(@plus, newE,C(1:end-1));
no_allcl = nClsOrig(end)+C(end-1);


