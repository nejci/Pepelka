function [clkd]=clust2kd(clust,n,m)
% [clkd]=clust2kd(clust,n,m)
%
% Convert clust vector as returned by clusot2d to
% a matrix decribing the clustering of neurons for
% maps given in kd format.
%
% Input:
% clust - (vector)
% n,m - (int) map dimensions (n x m)
%
% Output:
% clkd - (matrix) 
% 
% $Id: clust2kd.m 940 2005-11-04 17:57:33Z dome $
% D. Brugger, 25 July 2005
% util/clust2kd.m

clkd = [];
% number of neurons
num_n = n*m;
for k=1:n
  clkd = [clkd; clust(k:n:num_n)];
end
