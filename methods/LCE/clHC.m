function CR = clHC(S, K)
%==========================================================================
% FUNCTION: CR = clHC(S, K)
% DESCRIPTION: This function performs the final clustering using Hierarchical
%              algorithms (Single-Linkage:SL, Complete-Linkage:CL and
%              Average-Linkage:AL) as consensus methods.
%
% INPUTS:     S = an N-by-N similarity matrix
%             K = the prefered number of clusters
%
% OUTPUTS: CR = an N-by-3 matrix of clustering results from SL, CL and AL,
%               respectively
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

CR = [];
d = stod(S); %convert similarity matrix to distance vector
% single linkage
Z = linkage(d,'single');
CR = [CR cluster(Z,'maxclust',K)];
% complete linkage
Z = linkage(d,'complete');
CR = [CR cluster(Z,'maxclust',K)];
% average linkage
Z = linkage(d,'average');
CR = [CR cluster(Z,'maxclust',K)];