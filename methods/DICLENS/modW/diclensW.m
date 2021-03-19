function [labelsCons,Kcons,time,differ] = diclensW(labelsEns,weights,K)

% [labelsCons, Kcons, time, differ] = diclensW(labelsEns, weights, K)
%--------------------------------------------------------------------------
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number
% Consensus function for combining ensemble of clusterings.
%
% Based on Java implementation: http://www.cs.umb.edu/~smimarog/diclens/.
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns     (matrix)    clustering ensemble; each COLUMN corresponds
%							  to one clustering (standard in Pepelka).
%                 (string)    name of the text file containing cluster
%                             labels. One clustering per ROW or in a binary
%                             format (each row corresponds to one cluster)
%   weights       (vector)    weights for each partition if labelsEns is matrix; 
%                             vector of weights for each cluster if
%                             labelsEns is string;
%                             if empty, we will use a vector of ones.
%   K             (scalar)    optional - number of clusters in consensus
%                             partition is automatically determined by
%                             default. However, user can force this
%                             number by this parameter. In case it is not
%                             possible, ignore this attempt.
%--------------------------------------------------------------------------
% OUTPUTS:
%   labelsCons    (vector)    ensemble clustering result - labels of
%                             consensus clustering in a column vector
%	Kcons         (scalar)	number of clusters in consensus clustering
%	time			(scalar)	execution time in seconds
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Mimaroglu, S. & Aksehirli, E.
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number.
% Computational Biology and Bioinformatics, IEEE/ACM Transactions on,
% 2012, 9, 408-420
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 9-September-2013 by Nejc Ilc
%
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
    PI = str2num(file);
else
    PI = labelsEns';
end

tic;

% M = number of clusterings
% N = number of data samples
[M,N] = size(PI);

% save the positions of NaNs
nanmask = isnan(PI);

nClustL = zeros(M,1);
for iM = 1:M
    lbl = PI(iM,:);
    [u,~,iB] = unique(lbl);
    nClustL(iM) = length(u)-sum(nanmask(iM,:));
    
    % repair labels that are not in proper sequential form
    if ~isequal(u, 1:nClustL(iM))
        fprintf(1,'Labels in line %d will be repaired!\n',iM);
        PI(iM,:) = iB;
    end
end
PI(nanmask) = nan;

nClust = sum(nClustL);

% Construct binary representation of PI and save it to P.
% Each cluster is a row with ones on indeces of containing data points.
P = false(nClust,N);
sep=[0;cumsum(nClustL)];

for p=1:M
    Ptmp = false(nClustL(p),N);
    
    step = nClustL(p);
    linInd = PI(p,:) + (0:step:step*N-1);
    linInd(isnan(linInd))=[];
    Ptmp(linInd) = 1;
    P( (sep(p)+1) : sep(p+1) ,:) = Ptmp;
end

%--------------------------------------------------------------------------
% Construct graph
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Graph is represented by sparse matrix G, where each entry determine an
% edge and weight at the same time. Weights of the edges are all nonzero
% entries in the lower triangle of the N-by-N sparse matrix G.
% G(1,2) = 0.5 means that there is edge
% between node 1 and node 2 and is weighted by 0.5.
%--------------------------------------------------------------------------

% compute evidence accumulation (co-occurence) matrix - each entry
% represents similarity between two data points


%AM = accumat(P);
%AM = computeCO(P,M,'binary','full',[]).*M;
AM = computeCO(PI',M,'labels','full',weights).*M;

% For speed, precompute number of ones in rows of P
nElem_P = sum(P,2);

G = sparse(nClust,nClust);
for p=1:M-1
    cInd = sep(p+1)+1:nClust;
    gInd = sep(p)+1:sep(p+1);
    
    for Cik = gInd
        for Cjl = cInd
            % we switch indeces to produce lower triangular matrix G
            %G(Cjl,Cik) = computeECS2(P(Cik,:),P(Cjl,:),AM);
            nElem = nElem_P(Cik)*nElem_P(Cjl);
            if nElem >0
                sumECS = sum(sum(AM(P(Cik,:),P(Cjl,:))));
                sumECS = sumECS / nElem;
                G(Cjl,Cik) = sumECS;
            end
        end
    end
end

% costs of edges and weights have inverse relationship:
maxW = max(max(G)) + 1;

V = nchoose2(1:nClust); % using FEX #20144 by Jos
I = V(:,1)';
J = V(:,2)';
S = maxW - G(J+nClust.*(I-1));  % Compute a linear index with vector operations
S(S==maxW)=0; % preserve zeroes (absence of an edge in graph)
Gtrans = sparse(J,I,S,nClust,nClust);


[SMST] = graphminspantree(Gtrans, 'Method', 'Prim');
%view(biograph(SMST,[],'ShowArrows','off','ShowWeights','on'))

% sort edges in ascending order, from minimal to maximal cost
[Eind_i, Eind_j,E_val] = find(SMST); % num of edges in MST = num of nodes+1
E_ind =[Eind_i, Eind_j];
[E_val,pInd] = sort(E_val,'ascend');
E_ind = E_ind(pInd,:);

E_len=length(E_ind);

metaClusterGraph = sparse(nClust, nClust);
avgICS = zeros(1,nClust);
avgECS = zeros(1,nClust);
compHist = zeros(nClust,nClust); % history of component's labels

finalClusters_hist = cell(1,E_len);

for edgeInd = 1 : E_len
    
    % find connected components in metaClusterGraph
    [numComp, comp] = graphconncomp(metaClusterGraph, 'Directed', false);
    compHist(numComp,:) = comp;
    if numComp > nClust
        continue;
    end
    
    finalClusters = majorityVoting(comp, numComp, P);
    finalClusters_hist{numComp} = finalClusters;
    
    nFin = size(finalClusters,1);
    nElem_fC = sum(finalClusters,2);
    
    
    tmp=0;
    for cInd = 1:nFin
        %tmp = tmp + computeICS2(finalClusters(cInd,:), AM);
        cLen = nElem_fC(cInd);
        if cLen > 1
            npair = cLen * (cLen -1); % to save time, do not divide by 2, because in the next line, we will have to multiply by 2.
            tmp = tmp + sum(sum(AM(finalClusters(cInd,:),finalClusters(cInd,:))))/npair;
        end
    end
    avgICS(numComp) = tmp / nFin;
    
    tmp=0;
    for cInd1 = 1:nFin-1
        for cInd2 = cInd1+1:nFin
            %tmp = tmp + computeECS2(finalClusters(cInd1,:), finalClusters(cInd2,:), AM);
            nElem = nElem_fC(cInd1)*nElem_fC(cInd2);
            if nElem > 0
                sumECS = sum(sum(AM(finalClusters(cInd1,:),finalClusters(cInd2,:))));
                tmp = tmp + sumECS / nElem;
            end
        end
    end
    avgECS(numComp) = tmp / (nFin * (nFin-1) /2);
    
    metaClusterGraph(E_ind(edgeInd,1),E_ind(edgeInd,2)) = E_val(edgeInd);
end

% normalize ICS and ECS using min-max method;

avgECS = normalize_minmax(avgECS);
avgICS = normalize_minmax(avgICS);

% quality function is difference between normalized ICS and ECS
differ = avgICS - avgECS;

% Allow user to force the number of output clusters.
% If input parameter (optional) K is defined, try to find such cutpoint
% that will output desired number of clusters. It this is not possible,
% return the partition with the highest differ value.

% the solution is the partition of metacluster graph with the highest differ
% value
[~,cutPoint] = max(differ);

if exist('K','var') && ~isempty(K)
    [numClust,~] = cellfun(@size,finalClusters_hist);
    differ_K = differ;
    differ_K(K~=numClust) = -Inf;
    
    if all(isinf(differ_K))
       warning('pplk:diclens',['Cannot force the number of output clusters to K=',num2str(K),'. It has been automatically determined.']); 
    else
       [~,cutPoint] = max(differ_K);
    end
end

finalClusters = finalClusters_hist{cutPoint};
[~,labelsCons] = max(finalClusters,[],1);
    
time = toc;
labelsCons = labelsCons';
Kcons = size(finalClusters,1);



