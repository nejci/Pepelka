function V = cleval(X, CR, methods, truelabels)
%==========================================================================
% FUNCTION: V = cleval(X, CR, methods, truelabels)
% DESCRIPTION: This function computes validity scores for clustering
%              results and display a comparison bar chart
%
% INPUTS:   X = a dataset, rows of X correspond to observations; columns
%               correspond to variables (exclude class labels!!)
%          CR = a matrix of clustering results (exclude row and cloumn headers!!)
%     methods = a set of legend strings (clustering methods) to be shown in
%               the bar chart
%  truelabels =(optional) known cluster labels for each data points (N-by-1 vector) 
%
% OUTPUTS:  V = matrix of cluster validity scores
%
% NOTE1: format of 'V' ==> each column refers to each clustering method
%                          row2 = Compactness (CP),
%                          row3 = Davies-Bouldin Index (DB),
%                          row4 = Dunn Index,
%                          row5 = Adjust Rand Index (AR),
%                          row6 = Rand Index (RI),
%                          row7 = Classification Accuracy (CA)
%        The last three rows (AR, RI and CA) will be displayed only when a user specify 'truelabels' argument
% NOTE2: CP and DB: low values indicate good cluster structures
%        Dunn, AR, RI and CA: large values indicate better cluster quality
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

V = [];
for i=1:size(CR,2)
    cp = valid_compactness(X, CR(:,i)); %compute Compactness index
    [db, dunn] = valid_DbDunn(X,CR(:,i)); %compute Davies-Bouldin index and Dunn index
    if exist('truelabels') %compare result with known cluster labels
        [AR,RI]=valid_RandIndex(CR(:,i),truelabels); % compute Rand index and its variations
        ca = valid_CA(CR(:,i), truelabels); %compute Classification Accuracy
        V = [V [cp;db;dunn;AR;RI;ca]];
    else
        V = [V [cp;db;dunn]];
    end
end

if exist('truelabels')
    xlabel = {'CP';'DB';'Dunn';'AR';'RI';'CA'};
else
    xlabel = {'CP';'DB';'Dunn'};
end
showBar(V, xlabel, methods);

% format result
V = [['Validity Index' methods]; [xlabel num2cell(V)]];