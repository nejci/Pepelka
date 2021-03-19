function showBar(V, xlabel, legends)
%==========================================================================
% FUNCTION: showBar(V, xlabel, legends)
% DESCRIPTION: This function shows a bar chart for cluster validity comparison
%
% INPUTS: V = matrix of cluster validity scores (exclude row and column headers!!)
%    xlabel = labels for x-axis
%   legends = a set of legend strings
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

bar(V);
legend(legends,'Location','NorthEastOutside');
set(gca,'XTickLabel',xlabel)
set(gcf,'Name','Cluster Validity Comparison')