function [labels]=cluster_data(sM,td,clust)
% [labels]=cluster_data(sM,td,clust)
%
% cluster_data computes the predicted cluster labels
% for training data td using clustering clust for the
% som map units.
%
% Inputs:
% sM - (som_struct) trained som map
% td  - (som_struct) training data
% clust - (vector) clustering result, where clust(i)=k,
%         if map-neuron i belongs to cluster k
%
% Output:
% label - (vector) predicted cluster labels, e.g
%                  label(i)=k means that input vector i is
%                  assigned to cluster k. label(i)=0 means that
%                  the input vector was not assigned to a cluster.
% 
% $Id: cluster_data.m 673 2005-05-04 15:13:22Z dome $
% D. Brugger, 02 May 2005
% util/cluster_data.m

bmus = som_bmus(sM,td);
labels = clust(bmus)';
% convention: label 0 means 'unassigned'
labels = labels - 1;
