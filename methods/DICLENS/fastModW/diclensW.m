function [labelsCons, Kcons, time] = diclensW(labelsEns,weights,K)

% [labelsCons, Kcons, time] = DICLENSW(labelsEns,weights,K)
%--------------------------------------------------------------------------
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number
% Consensus function for combining ensemble of clusterings.
%
% Based on Java implementation: http://www.cs.umb.edu/~smimarog/diclens/.
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns		(matrix)	clustering ensemble; each COLUMN corresponds
%								to one clustering (standard in Pepelka).
%					(string)	name of the text file containing cluster
%								labels. One clustering per ROW or in a
%								binary format (each row corresponds to one
%								cluster)
%   weights       (vector)    weights for each partition if labelsEns is matrix; 
%                             vector of weights for each cluster if
%                             labelsEns is string;
%                             if empty, we will use a vector of ones.
%   K               (scalar)    optional - number of clusters in consensus
%                               partition is automatically determined by
%                               default. However, user can force this
%                               number by this parameter. In case it is not
%                               possible, ignore this attempt.
%--------------------------------------------------------------------------
% OUTPUTS:
%   labelsCons		(vector)	ensemble clustering result - labels of
%								consensus clustering in a column vector
%	Kcons			(scalar)	number of clusters in consensus clustering
%	time			(scalar)	execution time in seconds
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2016  Nejc Ilc
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Mimaroglu, S. & Aksehirli, E.
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number.
% Computational Biology and Bioinformatics, IEEE/ACM Transactions on,
% 2012, 9, 408-420
%
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 18-May-2016 by Nejc Ilc
%   - optimization for speed (mex)
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ischar(labelsEns)
    %--------------------------------------------------------------------------
    % PARSE
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Each column corresponds to one data sample.
    % Each row represents a clustering and each column shows cluster membership of
    % an object in the cluster. Order of the rows is not important.
    % Objects are separated with commas.
    % Unknown values are given as 'N'.
    %--------------------------------------------------------------------------
    file = fileread(labelsEns);
    file = strrep(file,'N','NaN');
    
    % PI is an representation of multiple clusterings in a matrix form.
    labelsEns = str2num(file);
    labelsEns = labelsEns';
end

ticID = tic();

% M = number of clusterings
% N = number of data samples
[N,M] = size(labelsEns);

% save the positions of NaNs
nanmask = isnan(labelsEns);

nClustL = zeros(M,1);
for iM = 1:M
    lbl = labelsEns(:,iM);
    [u,~,iB] = unique(lbl);
    labelsEns(:,iM) = iB;
    nClustL(iM) = length(u)-sum(nanmask(:,iM));
end
labelsEns(nanmask) = nan;


% Construct binary representation of PI and save it to P.
% Each cluster is a row with ones on indeces of containing data points.
[labelsEns, nClust, sep] = relabelCl(labelsEns);
[P,nElem_P] = computeP_mex(labelsEns, nClust);

% Compute co-occurence matrix AM
AM = computeCO(labelsEns,M,'labels','full',weights).*M;
AM(1:N+1:N^2) = M;

% Compute similarity graph Gtrans
% How much of memory can we take in advance for construction of graph
memoryLimitGB = 1;
[Gi,Gj,Gv] = computeG_mex(P,AM,sep,nElem_P,memoryLimitGB);
Gtrans = sparse(Gi,Gj,Gv,nClust,nClust,length(Gv));

% Compute Minimum Spanning Tree on graph
[Eind_i Eind_j E_val] = mst_mex(Gtrans,'prim','matrix',0);
clear Gtrans Gi Gj Gv;
% sort edges in ascending order, from minimal to maximal cost
[E_val,pInd] = sort(E_val,'ascend');
Eind_i = Eind_i(pInd);
Eind_j = Eind_j(pInd);
E_len = length(Eind_i);
% [Eind_i, Eind_j, E_val]

% Enter main loop
[avgICS, avgECS, numClustHist] = DiclensMainLoop_mex(P,nElem_P,AM,Eind_i,Eind_j,E_val);

% normalize ICS and ECS using min-max method;
avgECS = normalize_minmax(avgECS);
avgICS = normalize_minmax(avgICS);

% quality function is difference between normalized ICS and ECS
differ = avgICS - avgECS;

% Allow user to force the number of output clusters.
% If input parameter (optional) K is defined, try to find such cutpoint
% that will output desired number of clusters. It this is not possible,
% return the partition with the highest differ value.
% The solution is the partition of metacluster graph with the highest differ
% value
[~,cutPoint] = max(differ(2:end)); % skip first element (empty element)
cutPoint = cutPoint+1;

if exist('K','var') && ~isempty(K)
    
    differ_K = differ;
    differ_K(K~=numClustHist) = -Inf;
    
    if all(isinf(differ_K))
        %warning('pplk:diclens',['Cannot force the number of output clusters to K=',num2str(K),'. It has been automatically determined.']);
    else
        [~,cutPoint] = max(differ_K);
    end
end

% For space savings, we do not store all finalClusters. Instead we compute
% once again only the final solution, which maximizes differ.
edgeInd = E_len - cutPoint +2;
% Create meta-cluster graph with cutPoint edges
mG_i = [Eind_i(1:edgeInd-1); Eind_j(1:edgeInd-1)];
mG_j = [Eind_j(1:edgeInd-1); Eind_i(1:edgeInd-1)];
mG_v = [E_val(1:edgeInd-1); E_val(1:edgeInd-1)];
metaClusterGraph = sparse(mG_i, mG_j,mG_v,nClust,nClust,2*edgeInd-1);

% find connected components in metaClusterGraph
[comp, compSizes] = components_mex(metaClusterGraph);

% Vote for meta clusters.
[finalClusters,~] = majorityVoting_mex(comp,compSizes,P,nElem_P,N);

% Find final clustering partition
[~,labelsCons] = max(finalClusters,[],1);

time = toc(ticID);

labelsCons = labelsCons';
Kcons = numClustHist(cutPoint);%size(finalClusters,1);



