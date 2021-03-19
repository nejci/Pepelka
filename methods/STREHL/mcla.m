% function cl = mcla(cls,k)
%
% DESCRIPTION
%  Performs MCLA for CLUSTER ENSEMBLES
%
% Copyright (c) 1998-2002 by Alexander Strehl

function cl = mcla(cls,k,weights)

if ~exist('k','var') || isempty(k)
    k = max(max(cls));
end

clb = clstoclbs(cls);
cl_lab = clcgraph(clb,k,'simbjac');
for i=1:max(cl_lab),
   matched_clusters = cl_lab==i;
   clb_cum(i,:) = mean(clb(matched_clusters,:),1);
end;
cl = clbtocl(clb_cum);

