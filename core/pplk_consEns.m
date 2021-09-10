function [labelsCons,numClust,moreInfo] = pplk_consEns(labelsEns,K,method,params)
% [labelsCons,numClust,moreInfo] = pplk_consEns(labelsEns,K,method,params)
% Computes consensus partition from partitions in the ensemble.
%
% INPUTS
%   labelsEns       
%       Cluster ensemble - labels are stored in the columns of an N-by-E
%       matrix, where N is number of data points and E is number of
%       clusterings.
% 
%   K
%       Number of final clusters; can be empty for some methods, which can
%       automatically determine K.
% 
%   method
%       Consensus function name, one of the following:
% 
%           'BCE'
%               P. Wang, C. Domeniconi, and K. Laskey, "Nonparametric
%               Bayesian Clustering Ensembles," in Machine Learning and
%               Knowledge Discovery in Databases, vol. 6323, J. Balcázar,
%               F. Bonchi, A. Gionis, and M. Sebag, Eds. Springer Berlin /
%               Heidelberg, 2010, pp. 435-450.
%           'COMUSA'
%               S. Mimaroglu and E. Erdil, "Combining multiple clusterings
%               using similarity graph," Pattern Recognition, vol. 44, no.
%               3, pp. 694-703, Mar. 2011.
%           'DICLENS'-[''|'W'|'ORIG']
%               S. Mimaroglu and E. Aksehirli, "DICLENS: Divisive
%               Clustering Ensemble with Automatic Cluster Number,"
%               IEEE/ACM transactions on computational biology and
%               bioinformatics, vol. 9, no. 2, pp. 408-420, Sep. 2012.
%           'EAC'-['SL'|'CL'|'AL'|'WL']
%               A. L. N. Fred and A. K. Jain, "Combining Multiple
%               Clusterings Using Evidence Accumulation," IEEE Transactions
%               on Pattern Analysis and Machine Intelligence, vol. 27, no.
%               6, pp. 835-850, 2005.
%           'HBGF'
%               X. Z. Fern and C. E. Brodley, "Solving cluster ensemble
%               problems by bipartite graph partitioning," in Proceedings
%               of the twenty-first international conference on Machine
%               learning, 2004, p. 36.
%           'HUANG'-['SL'|'CL'|'AL'|'GPMGLA']
%               Dong Huang, Jian-Huang Lai, Chang-Dong Wang. Combining
%               Multiple Clusterings via Crowd Agreement Estimation and
%               Multi- Granularity Link Analysis. Neurocomputing, 2014.
%           'LCE'-['CTS'|'SRS'|'ASRS']-['SL'|'CL'|'AL']
%               N. Iam-On, T. Boongoen, and S. Garrett, "LCE: a link-based
%               cluster ensemble method for improved gene expression data
%               analysis," Bioinformatics, vol. 26, no. 12, pp. 1513-1519,
%               2010.
%           'PAC'-['SL'|'CL'|'AL'|'WL']
%               X. Wang, C. Yang, and J. Zhou, "Clustering aggregation by
%               probability accumulation," Pattern Recognition, vol. 42,
%               pp. 668-675, 2009.
%           'STREHL'-['CSPA'|'HGPA'|'MCLA']['-W']
%               A. Strehl and J. Ghosh, "Cluster ensembles - a knowledge
%               reuse framework for combining multiple partitions," The
%               Journal of Machine Learning Research, vol. 3, pp. 583-617,
%               2003.
%           'WEA'-['SL'|'CL'|'AL'|'WL']
%               S. Vega-Pons, J. Ruiz-Shulcloper, and A. Guerra-Gandón,
%               "Weighted association based methods for the combination of
%               heterogeneous partitions," Pattern Recognition Letters,
%               vol. 32, no. 16, pp. 2163-2170, Dec. 2011.
%           'WEAC'-['SL'|'CL'|'AL'|'WL']
%               F. J. Duarte, A. L.N.Fred, A. Lourenco, and M. F.
%               Rodrigues, "Weighting Cluster Ensembles in Evidence
%               Accumulation Clustering," in 2005 Portuguese Conference on
%               Artificial Intelligence, 2005, vol. 00, pp. 159-167.
%           'ZHOU'-['VOTE'|'WVOTE'|'SVOTE'|'SWVOTE']
%               IMPORTANT! Number of clusters in ensemble have to equal
%               number in consensus. Z.-H. Zhou and W. Tang, "Clusterer
%               ensemble," Knowledge-Based Systems, vol. 19, no. 1, pp.
%               77-83, 2006.
%
%   params
%       Parameters structure; can be empty or non-existent for defaults.
%       Used fields:
%
%       BCE
%           BCE_Ki        
%              Number of clusters in ensemble.
%           BCE_alphaInit 
%             Init values for alpha.
%           BCE_betaInit  
%             Init values for beta.
%           BCE_lapParam  
%               Laplacian smoothing parameter.
%       COMUSA
%           COMUSA_relaxation 
%               Relaxation of maximum edge weight constraint. 0.5 means 50%
%               relaxation and thus bigger final clusters [0].
%       HUANG
%           HUANG_betaOrINCAI
%           HUANG_alpha
%       LCE
%           LCE_dc        
%               Decay factor [0.8].
%           LCE_R         
%               Number of iters for SimRank algorithm [5].
%       PAC
%           PAC_dim       
%               Dimensionality of the original data.
%       WEA
%           WEA_data      
%               Original data or (dis)similarity matrix.
%           WEA_dataMode  
%               'data', 'dist' or 'sim'.
%           WEA_dataDist  
%               If dataMode is 'data', provide with distance type between 
%               data.
%           WEA_normalize 
%               Whether to use normalized weights [1].
%       WEAC and DICLENS-W
%           WEAC_data
%           WEAC_CVI
%           WEAC_CVImat
%           WEAC_unifyMeth
%           WEAC_reduceMeth
%           WEAC_reduceDim
%           WEAC_weightMeth
%           WEAC_weightMode
%           WEAC_weights
%
% 
% OUTPUTS
%   labelsCons
%       Ensemble clustering result - labels.
%
%   numClust
%       Number of clusters in the consensus partition.
%
%   TODO Moreinfo explenation
% 
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka


callDir=chdir(pplk_homeDir());

moreInfo = [];

if ~exist('params','var') || isempty(params)
    params=pplk_setParamsDefault();
end

M=size(labelsEns,2);

% Check labelsEns for consistency (if the labels are in form
% 1:length(unique(labels_i))) and fix if needed
ensMember_K = zeros(1,M);
for iM = 1:M
    lbl = labelsEns(:,iM);
    [u,~,iB] = unique(lbl);
    ensMember_K(iM) = length(u);
    % repair labels that are not in proper sequential form
    labelsEns(:,iM) = iB;
end

% split the string that describes methods for composing similarity
% matrices and final clustering (see LCE)
method = upper(method);
methodsTok=regexp(method,'-','split');
nTok=length(methodsTok);

%--------------------------------------------------------------------------
% PRM step
% determine if weigted mode is ON
weightedModeOn = 0;
if strcmpi(methodsTok{1},'WEAC')
    weightedModeOn = 1;
elseif strcmpi(methodsTok{1},'DICLENS')
    if nTok > 1 && any(strcmpi(methodsTok{2},{'W','WALT'}))
        weightedModeOn = 1;
    end
elseif strcmpi(methodsTok{1},'STREHL')
    if nTok > 2 && strcmpi(methodsTok{3},{'W'})
        weightedModeOn = 1;
    end
end

% Weighted mode - PRM step, compute weights
if weightedModeOn
    
    % some defaults
    weightMeth = 'wMean2';
    weightMode = '';
    data = [];
    
    if ~isfield(params,'WEAC_weights')
        weights = [];
    else
        weights = params.WEAC_weights;
    end
    
    % If weights are not known, compute them
    if isempty(weights)
        PRMopt = [];
        
        if ~isfield(params,'WEAC_CVI')
            error('List of cluster validity indices params.WEAC_CVI is missing.');
        else
            CVI = params.WEAC_CVI;
        end
        
        if ~isfield(params,'WEAC_CVImat')
            CVImat = [];
        else
            CVImat = params.WEAC_CVImat;
        end
        
        if ~isempty(CVImat)
            PRMopt.CVImat = CVImat;
        else
            if ~isfield(params,'WEAC_data')
                error('Method WEAC requires params.WEAC_data to compute cluster validity indices.');
            else
                data = params.WEAC_data;
            end
        end
        
        if ~isfield(params,'WEAC_unifyMeth')
            unifyMeth = [];
        else
            unifyMeth = params.WEAC_unifyMeth;
        end
        
        if ~isfield(params,'WEAC_reduceMeth')
            reduceMeth = [];
        else
            reduceMeth = params.WEAC_reduceMeth;
        end
        
        if ~isfield(params,'WEAC_reduceDim')
            reduceDim = [];
        else
            reduceDim = params.WEAC_reduceDim;
        end
        PRMopt.reduceDim = reduceDim;
        
        if isfield(params,'WEAC_weightMeth')
            weightMeth = params.WEAC_weightMeth;
        end
        
        if isfield(params,'WEAC_weightMode')
            weightMode = params.WEAC_weightMode;
        end
        
        [weights,~,featsInd] = ...
            pplk_partitionRelevance(data,labelsEns,CVI,...
            unifyMeth,reduceMeth,weightMeth,weightMode,PRMopt);
        moreInfo.featsInd = featsInd;
    end
    % If all weights elements are zero, set them to 1.
    % All-zeros may cause error in weighted consensus function.
    if all(weights==0)
        weights = ones(size(weights));
        warning('All weights are zero! Setting them to 1.');
    end
end



switch methodsTok{1}
    
    case 'BCE'
        %----------------------------------------------------------------------
        % BCE - Bayesian Cluster Ensembles
        %        Wang, Shan & Banerjee, 2011
        %----------------------------------------------------------------------
        Ki = ensMember_K;
        alphaInit = [];
        betaInit = [];
        lapParam = [];
        if isfield(params,'BCE_Ki')
            Ki = params.BCE_Ki;
        end
        if isfield(params,'BCE_alphaInit')
            alphaInit = params.BCE_alphaInit;
        end
        if isfield(params,'BCE_betaInit')
            betaInit = params.BCE_betaInit;
        end
        if isfield(params,'BCE_lapParam')
            lapParam = params.BCE_lapParam;
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'BCE',filesep,'fast']);
        
        [labelsCons,numClust] = BCE(labelsEns,K,Ki,alphaInit,betaInit,lapParam);
        
        chdir(oldPath);
        
    case 'BCE_ORIG'
        %----------------------------------------------------------------------
        % BCE - Bayesian Cluster Ensembles
        %        Wang, Shan & Banerjee, 2011
        %----------------------------------------------------------------------
        Ki = ensMember_K;
        alphaInit = [];
        betaInit = [];
        lapParam = [];
        if isfield(params,'BCE_Ki')
            Ki = params.BCE_Ki;
        end
        if isfield(params,'BCE_alphaInit')
            alphaInit = params.BCE_alphaInit;
        end
        if isfield(params,'BCE_betaInit')
            betaInit = params.BCE_betaInit;
        end
        if isfield(params,'BCE_lapParam')
            lapParam = params.BCE_lapParam;
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'BCE',filesep,'orig']);
        
        [labelsCons,numClust] = BCE(labelsEns,K,Ki,alphaInit,betaInit,lapParam);
        
        chdir(oldPath);
        
    case 'EAC'
        %----------------------------------------------------------------------
        % EAC - Evidence Accumulation Clustering
        %        Fred & Jain, 2005
        %----------------------------------------------------------------------
        oldPath=chdir(['..',filesep,'methods',filesep,'EAC']);
        if nTok > 1
            consFun = methodsTok(2:end);
        else
            consFun = [];
        end
        
        [labelsCons,numClust] = EAC(labelsEns,K,consFun);
        
        chdir(oldPath);
        
        
    case 'HUANG'
        %----------------------------------------------------------------------
        % HUANG - Crowd Agreement Estimation and Multi-Granularity Link Analysis
        %        Huang et al., 2014
        %----------------------------------------------------------------------
        if ~isfield(params,'HUANG_betaOrINCAI')
            params.HUANG_betaOrINCAI = [];
        end
        if ~isfield(params,'HUANG_alpha')
            params.HUANG_alpha = [];
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'HUANG']);
        
        if nTok > 1
            consFun = methodsTok(2:end);
        else
            consFun = 'GPMGLA';
        end
        
        if strcmpi(consFun,'GPMGLA')
            [labelsCons, numClust] = HUANG_GPMGLA(labelsEns, K, ...
                params.HUANG_betaOrINCAI, params.HUANG_alpha);
        else
            [labelsCons, numClust] = HUANG_WEAC(labelsEns, K, consFun,...
                params.HUANG_betaOrINCAI);
        end
        chdir(oldPath);
        
    case 'WEAC'
        %----------------------------------------------------------------------
        % WEAC - Weighted Evidence Accumulation Clustering
        %        Duarte et al., 2005
        %        Pepelka extensions, 2013
        %----------------------------------------------------------------------
        
        % TODO: test adaptation for multi consFun
        
        oldPath=chdir(['..',filesep,'methods',filesep,'EAC']);
        
        if nTok > 1
            consFun = methodsTok(2:end);
        else
            consFun = [];
        end
        
        WEACopt = [];
        WEACopt.weights = weights;
        
        % note: data, CVI and weightMeth are not used in this call, because
        % weights are known.
        [labelsCons,numClust] = WEAC([],labelsEns,K,consFun,[],[],WEACopt);
        
        chdir(oldPath);
        
    case 'PAC'
        %----------------------------------------------------------------------
        % PAC - Probability Accumulation Clustering
        %        Wang et al., 2009
        %----------------------------------------------------------------------
        if ~isfield(params,'PAC_dim')
            error('Method PAC requires dimensionality of data to be known.');
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'PAC']);
        
        if nTok > 1
            consFun = methodsTok(2:end);
        else
            consFun = [];
        end
        
        [labelsCons,numClust] = PAC(labelsEns,params.PAC_dim,K,consFun);
        
        chdir(oldPath);
        
    case 'WEA'
        %----------------------------------------------------------------------
        % WEA - Weighted Evidence Accumulation
        %        Vega-Pons et al., 2009, 2011
        %----------------------------------------------------------------------
        
        if ~isfield(params,'WEA_data')
            error('Method WEA requires data as parameter.');
        end
        if ~isfield(params,'WEA_dataMode')
            params.WEA_dataMode = [];
        end
        if ~isfield(params,'WEA_dataDist')
            params.WEA_dataDist = [];
        end
        if ~isfield(params,'WEA_normalize')
            params.WEA_normalize = [];
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'WEA']);
        
        if nTok > 1
            consFun = methodsTok(2:end);
        else
            consFun = [];
        end
        [labelsCons,numClust] = WEA(labelsEns,params.WEA_data, ...
            params.WEA_dataMode,params.WEA_dataDist,K,consFun,...
            params.WEA_normalize);
        
        chdir(oldPath);
        
    case 'STREHL'
        %------------------------------------------------------------------
        % CSPA, HGPA and MCLA
        %        Strehl & Ghosh, 2002
        %------------------------------------------------------------------
        %Strehl's CSPA, HGPA and MCLA
        % we choose result with the highest NMI
        
        if nTok < 2
            methodsTok{2}='CSPA';
        end
        
        
        if nTok > 1
            consFun = methodsTok{2};
        else
            consFun = [];
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'STREHL']);
        
        if weightedModeOn
            labelsCons = clusterensemble(labelsEns',K,consFun,weights)';
        else
            labelsCons = clusterensemble(labelsEns',K,consFun)';
        end
        
        numClust = length(unique(labelsCons));
        chdir(oldPath);
        
    case 'LCE_ORIG'
        %------------------------------------------------------------------
        % LCE  (Link-Based Cluster Ensemble)
        %       Iam-On & Garrett, 2010
        %------------------------------------------------------------------
        if nTok < 3
            methodsTok{3}='SL';
        end
        if nTok < 2
            methodsTok{2}='CTS';
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'LCE',filesep,'orig']);
        labelsCons = runLCE(labelsEns,K,methodsTok{2},methodsTok{3},params.LCE_dc,params.LCE_R);
        numClust = length(unique(labelsCons));
        chdir(oldPath);
        
    case 'LCE'
        %------------------------------------------------------------------
        % LCE  (Link-Based Cluster Ensemble)
        %       Iam-On & Garrett, 2010
        %      Optimized for speed - Pepelka.
        %------------------------------------------------------------------
        if nTok < 3
            consFun = 'SL';
        else
            consFun = methodsTok(3:end);
        end
        if nTok < 2
            refinedMeth = 'CTS';
        else
            refinedMeth = methodsTok{2};
        end
        
        oldPath=chdir(['..',filesep,'methods',filesep,'LCE']);
        labelsCons = runLCE(labelsEns,K,refinedMeth,consFun,params.LCE_dc,params.LCE_R);
        numClust = length(unique(labelsCons));
        chdir(oldPath);
        
    case 'HBGF'
        %------------------------------------------------------------------
        % HBGF (Hybrid Bipartite Graph Formulation)
        %       Fern & Brodley, 2004
        %------------------------------------------------------------------
        oldPath=chdir(['..',filesep,'methods',filesep,'HBGF']);
        labelsCons = HBGF(labelsEns,K);
        numClust = length(unique(labelsCons));
        chdir(oldPath);
        
        
    case 'COMUSA'
        %------------------------------------------------------------------
        % COMUSA  (Combining multiple clusterings using similarity graph)
        %       Mimaroglu & Erdil, 2011
        %------------------------------------------------------------------
        oldPath=chdir(['..',filesep,'methods',filesep,'COMUSA']);
        [labelsCons,numClust] = comusa(labelsEns,params.COMUSA_relaxation);
        chdir(oldPath);
        
    case 'DICLENS'
        %------------------------------------------------------------------
        % DICLENS (Divisive Clustering Ensemble)
        %       Mimaroglu & Aksehirli, 2012
        %------------------------------------------------------------------
        
        % some experimental modes
        % 'DICLENS-ORIG': original SLOW version with bugs (Java port)
        % 'DICLENS'     : faster version without known bugs
        % 'DICLENS-W'   : weighted version based on fast DICLENS
        % 'DICLENS-WALT'   : weighted version based on alternative impl.
        
        if nTok < 2
            methodsTok{2}='';
        end
        
        switch(upper(methodsTok{2}))
            case ''
                oldPath=chdir(['..',filesep,'methods',filesep,'DICLENS',filesep,'fast']);
                [labelsCons,numClust] = diclens(labelsEns,K);
                
            case 'ORIG'
                oldPath=chdir(['..',filesep,'methods',filesep,'DICLENS',filesep,'orig']);
                [labelsCons, numClust] = diclens(labelsEns);
                
            case 'W'
                oldPath=chdir(['..',filesep,'methods',filesep,'DICLENS',filesep,'fastModW']);
                
                [labelsCons, numClust] = diclensW(labelsEns,weights,K);
                                
            case 'WALT'
                oldPath=chdir(['..',filesep,'methods',filesep,'DICLENS',filesep,'modW']);
                                
                [labelsCons, numClust] = diclensW(labelsEns,weights,K);
                
            otherwise
                error('Wrong mode for DICLENS consensus function.');
        end
        
        chdir(oldPath);
        
        
    case 'ZHOU'
        %----------------------------------------------------------------------
        % ZHOU (selective/weighted clusterer ensemble)
        %       Zhou & Tang, 2006
        %----------------------------------------------------------------------
        oldPath=chdir(['..',filesep,'methods',filesep,'ZHOU']);
        
        if nTok < 2
            [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',1,'weighted',1);
        else
            switch(methodsTok{2})
                case 'VOTE' % voting
                    [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',0,'weighted',0);
                case 'WVOTE' % weighted (NMI) voting
                    [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',0,'weighted',1);
                case 'SVOTE' % selective voting
                    [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',1,'weighted',0);
                case {'','SWVOTE'} % selective weighted voting
                    [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',1,'weighted',1);
                otherwise
                    warning('Pepelka:wrongMethod',['Method ', method,' does not exist! Using ZHOU-SWVOTE.']);
                    [labelsCons, numClust] = clusterer_ensemble_mod(labelsEns,K,'selective',1,'weighted',1);
            end
        end
        chdir(oldPath);
        
    otherwise
        error('Wrong consensus method name!');
end

% Check the number of clusters and output a warning if desired number does not
% match the number in consensus clustering.
if ~isempty(K) && any(numClust ~= K)
    warning('Pepelka:numberOfClusters',['Number of clusters (',num2str(numClust),') different from K (',num2str(K),').']);
    %fprintf(1,['pplk_consEns: number of clusters (',num2str(numClust),') different from K (',num2str(K),').']);
end

chdir(callDir);