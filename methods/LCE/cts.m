function S = cts(E, dc) 
%==========================================================================
% FUNCTION: S = cts(E, dc) 
% DESCRIPTION: A funciton for computing Connected-Triple Based Similarity matrix
%
% INPUT:  E = matrix of cluster ensemble
%        dc = decay factor, ranges [0,1]
%
% OUTPUT: S = CTS matrix
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

[E, nCls, C] = relabelCl(E); % re-labelling clusters in the ensemble E
wcl= weightCl(E,nCls);
S = cts_S_mex(E,wcl,C,dc);



