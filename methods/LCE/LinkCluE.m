function [CR,V] = LinkCluE(X, M, k, scheme, K, dcCTS, dcSRS, R, dcASRS, truelabels)
%==========================================================================
% FUNCTION: [CR,V] = LinkCluE(X, M, k, scheme, K, dcCTS, dcSRS, R, dcASRS, truelabels)
% DESCRIPTION: This function performs link-based cluster ensemble for the dataset X
%
% INPUTS:   X = a dataset, rows of X correspond to observations; columns
%               correspond to variables (exclude class labels!!)
%           M = the prefered number of base clusterings in the ensemble 
%           k = the prefered number of clusters in the base clusterings (should be greater than K)
%      scheme = cluster ensemble generating scheme (1 = Fixed k, 2 = Random k)
%           K = the number of clusters using in consensus functions (expected number of clusters)
%       dcCTS = decay factor (ranges [0,1]) for CTS method
%       dcSRS = decay factor (ranges [0,1]) for SRS method
%           R = the number of iterations for SimRank algorithm (don't need
%               for CTS)
%       dcASRS = decay factor (ranges [0,1]) for ASRS method
%  truelabels =(optional) known cluster labels for each data points (N-by-1 vector) 
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



% Argument checking ---------------
if nargin<9; dcASRS = 0.8; end;  
if nargin<8; R = 5; end;   
if nargin<7; dcSRS = 0.8; end; 
if nargin<6; dcCTS = 0.8; end;

row = size(X,1);
if row > 1000
    warning('The number of data points > 1000, it may take long time to process.');
end

mNames = {'cts', 'srs', 'asrs'};
dcs = {dcCTS, dcSRS, dcASRS};

if ~ismember(scheme, [1 2])
    error('"scheme" should be 1 or 2 for "Fixed k" and "Random k", respectively.')
end

for i = 1:length(dcs)
    dc = dcs{i};
    if (dc > 1) || (dc < 0)
        error('"dc" must be in the range of 0 and 1.')
    end
end

intVals = {M, k, K, R};

for i = 1:length(intVals)
    intVal = intVals{i};
    if (ceil(intVal) ~= floor(intVal)) || (intVal < 0)
        error([num2str(intVal) ' must be a positive integer.'])
    end
end
%-------------------------------


% Generate cluster ensemble
E = crEnsemble(X, M, k, scheme);


% Generating similarity matrices and performing consensus functions
CR = [];
for i = 1:length(mNames)
    mName = mNames{i};
    dc = dcs{i};
    disp(['Generating ' upper(mName) ' matrix...']);
    if ~strcmp(mName,'srs')
        S = feval(mName, E, dc); % CTS and ASRS
    else
        S = feval(mName, E, dc, R); % SRS
    end
    
    CR = [CR clHC(S, K)]; % perform consensus functions
end


% Evaluate quality of clutering results
methods = {'CTS-SL','CTS-CL','CTS-AL','SRS-SL','SRS-CL','SRS-AL','ASRS-SL','ASRS-CL','ASRS-AL'};
if ~exist('truelabels','var')
    V = cleval(X, CR, methods);
else
    V = cleval(X, CR, methods, truelabels);
end

% format result
CR = [methods; [num2cell(CR)]];