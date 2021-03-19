function [CR,V] = LinkCluETest
%==========================================================================
% FUNCTION: [CR,V] = LinkCluETest
% DESCRIPTION: A function for testing link-based cluster ensemble algorithms
% 
% OUTPUTS: CR = matrix of clustering results
%           V = matrix of cluster validity scores
%
% NOTE1: format of 'CR' ==> each column refers to each clustering method 
%                           e.g., 'CTS-SL' refers to CTS matrix with
%                           Single-Linkage algorithm.
%                           each row represents cluster labels for each data point
%        format of 'V'  ==> each column refers to each clustering method
%                           row2 = Compactness (CP),
%                           row3 = Davies-Bouldin Index (DB),
%                           row4 = Dunn Index,
%                           row5 = Adjust Rand Index (AR),
%                           row6 = Rand Index (RI),
%                           row7 = Classification Accuracy (CA)
%        The last three rows (AR, RI and CA) will be displayed only when a user specify 'truelabels' argument
% NOTE2: CP and DB: low values indicate good cluster structures
%        Dunn, AR, RI and CA: large values indicate better cluster quality
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

% identify all input arguments

%%%% for Four-Gaussian dataset %%%%%
X = load ('SampleData\FGD.csv'); %import Four-Gaussian data
K = 4; % the number of clusters in the final clustering (using in consensus functions)
truelabels = load ('SampleData\FGT.csv'); %import Four-Gaussian truelabels

%%%% for Leukemia dataset %%%%%
% X = load ('SampleData\LD.csv'); %import Leukemia data
% K = 2; % the number of clusters in the final clustering (using in consensus functions)
% truelabels = load ('SampleData\LT.csv'); %import Leukemia truelabels

M = 10; % the number of clusterings in ensemble

k = ceil(sqrt(size(X,1))); % the number of clusters in base clusterings

scheme = 1; % ensemble generating scheme (1 = fixed k, 2 = random k)

dcCTS = 0.8; % the decay factor fot CTS method

dcSRS = 0.8; % the decay factor fot SRS method

R = 5; % the number of iterations for SRS method (SimRank algorithm)

dcASRS = 0.8; % the decay factor fot ASRS method

% perform link-based cluster ensemble algorithm
[CR,V] = LinkCluE(X, M, k, scheme, K, dcCTS, dcSRS, R, dcASRS, truelabels); %truelabels is optional

% so, the function can be called:
% [CR,V] = LinkCluE(X, M, k, scheme, K, dcCTS, dcSRS, R); 