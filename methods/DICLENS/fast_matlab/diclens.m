function [labelsCons, Kcons, time] = diclens(labelsEns,K)

% [labelsCons, Kcons, time] = DICLENS(labelsEns,K)
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
% Version: 1.2
% Last modified: 15-July-2014 by Nejc Ilc
%   - optimization for speed (code from functions ECS and ICS is now in
%     main program)
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
    PI(iM,:) = iB;
    nClustL(iM) = length(u)-sum(nanmask(iM,:));
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
AM = computeCO(PI',M,'labels','full',[]).*M;
AM0 = AM; % for computing ICS
AM(1:N+1:N^2) = M; %sum(P,1); % !!! BUG??? similarity between the same point is what? M? Or should we count it (in case when point is not labeled in every partition -> toy example in the paper)

% For speed, precompute number of ones in rows of P
nElem_P = sum(P,2);

Gi = zeros(1,nClust);
Gj = zeros(1,nClust);
Gv = zeros(1,nClust);

Gidx = 1;
for p=1:M-1
    cInd = sep(p+1)+1:nClust;
    gInd = sep(p)+1:sep(p+1);
    
    for Cik = gInd
        for Cjl = cInd
            % we switch indeces to produce lower triangular matrix G
            nElem = nElem_P(Cik)*nElem_P(Cjl);
            if nElem >0
                sumECS = sum(sum(AM(P(Cik,:),P(Cjl,:)))) / nElem;
                Gi(Gidx) = Cik;
                Gj(Gidx) = Cjl;
                Gv(Gidx) = sumECS;
                Gidx = Gidx+1;
            end
        end
        
    end
end
Gi(Gidx:end) = [];
Gj(Gidx:end) = [];
Gv(Gidx:end) = [];
G = sparse(Gj,Gi,Gv,nClust,nClust,Gidx-1);

% costs of edges and weights have inverse relationship:
maxW = max(max(G)) + 1;
funSim = @(x) (maxW - x);
Gtrans = spfun(funSim,G);

% Using Graph BGL library
[Eind_i Eind_j E_val] = mst_mex(Gtrans'+Gtrans,'prim','matrix',0);

% sort edges in ascending order, from minimal to maximal cost
[E_val,pInd] = sort(E_val,'ascend');
Eind_i = Eind_i(pInd);
Eind_j = Eind_j(pInd);
E_len = length(Eind_i);
% [Eind_i, Eind_j, E_val]

avgICS = zeros(1,nClust);
avgECS = zeros(1,nClust);
finalClusters_hist = cell(1,nClust);


for edgeInd = 0:E_len-1
    
    % Add an edge from SMST into meta-cluster graph. Make it symmetric.
    mG_i = [Eind_i(1:edgeInd); Eind_j(1:edgeInd)];
    mG_j = [Eind_j(1:edgeInd); Eind_i(1:edgeInd)];
    mG_v = [E_val(1:edgeInd); E_val(1:edgeInd)];
    metaClusterGraph = sparse(mG_i, mG_j,mG_v,nClust,nClust,2*edgeInd);
    
    % find connected components in metaClusterGraph
    [comp, compSizes] = components_mex(metaClusterGraph);
    comp = comp';
    numComp = length(compSizes);
    
    % Vote for meta clusters.
    %finalClusters = majorityVoting(comp, numComp, P);
    
    % calculate votes
    votes = zeros(numComp, N);
    for compInd = 1:numComp
        metaCluster = (comp==compInd);
        votes(compInd,:)=sum(P(metaCluster,:),1);
    end
    
    % It depends on the order of rows in matrix votes when searching for max
    % value. If there is a tie, we choose one at random!
    randP = randperm(numComp);
    votes = votes(randP,:); % to simulate behaviour of Java implementation
    
    % assign data samples to the meta cluster with the highest vote rate
    finalClusters = false(numComp,N);
    [~, maxVotesInd] = max(votes,[],1);
    finalClusters( sub2ind([numComp,N],maxVotesInd,1:N) ) = 1;
    
    % remove empty clusters
    count = sum(finalClusters,2);
    finalClusters(count==0,:)=[];
    
    finalClusters_hist{numComp} = finalClusters;
    
    % Compute ICS and ECS of finalCLusters
    nElem_fC = sum(finalClusters,2);
    nFin = size(finalClusters,1);    
    ICS=0;
    for cInd = 1:nFin
        cLen = nElem_fC(cInd);
        if cLen > 1
            npair = cLen * (cLen -1); % to save time, do not divide by 2, because in the next line, we will have to multiply by 2.
            ICS = ICS + sum(sum(AM0(finalClusters(cInd,:),finalClusters(cInd,:))))/npair;
        end
    end
    ICS = ICS / nFin;
    
    ECS=0;
    for cInd1 = 1:nFin-1
        for cInd2 = cInd1+1:nFin
            nElem = nElem_fC(cInd1)*nElem_fC(cInd2);
            if nElem > 0
                sumECS = sum(sum(AM(finalClusters(cInd1,:),finalClusters(cInd2,:))));
                sumECS = sumECS / nElem;
                ECS = ECS + sumECS;
            end
        end
    end
    ECS = ECS / (nFin * (nFin-1) /2);
    
    avgICS(numComp) = ICS;
    avgECS(numComp) = ECS;
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
% The solution is the partition of metacluster graph with the highest differ
% value
[~,cutPoint] = max(differ(2:end)); % skip first element (empty element)
cutPoint = cutPoint+1;

if exist('K','var') && ~isempty(K)
    [numClust,~] = cellfun(@size,finalClusters_hist);
    differ_K = differ;
    differ_K(K~=numClust) = -Inf;
    
    if all(isinf(differ_K))
        %warning('pplk:diclens',['Cannot force the number of output clusters to K=',num2str(K),'. It has been automatically determined.']);
    else
        [~,cutPoint] = max(differ_K);
    end
end

% ??? Ce zelis privarcevati prostor, potem si shranjuj comp in tukaj klici
% finalClusters = majorityVoting(compHist(cutPoint,:), cutPoint, P);
if (isempty(finalClusters_hist{cutPoint}))
    fprintf(1,'!!!\n');
end
finalClusters = finalClusters_hist{cutPoint};
[~,labelsCons] = max(finalClusters,[],1);

time = toc;
labelsCons = labelsCons';
Kcons = size(finalClusters,1);



