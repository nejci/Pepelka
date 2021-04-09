function [validInt, list] = pplk_validInt(data, labels, methods, options)
% PPLK_VALIDINT 
% [validInt, list] = PPLK_VALIDINT(data, labels, methods, options) 
% Internal validity indices for estimating clustering quality and
% number of clusters.
%--------------------------------------------------------------------------
% INPUTS
%   data		: data that were clustered [nPatterns x nDimensions]
%
%   labels      : labels of data, result of a clustering task
%                 can be vector [nPatterns x 1] or matrix of labels [nPatterns x nNCs]
%
%	methods     : a cell containing one or more method's identifier.
%				  You can use multiple methods, i.e., {'DN','SIL','WR'}.
%				  If methods=[], all of them are computed.
%
%                 Note: indices that are designed ONLY for estimating number
%                 of clusters are indicated with *
%
% ID        Method name                     Optimum	Reference
% -------------------------------------------------------------------------
% 'APN'     Avg. proportion of non-overlap	min     Datta & Datta, 2003
% 'AD'      Average distance                min     Datta & Datta, 2003
% 'ADM'     Average distance between means  min     Datta & Datta, 2003
% 'BHI'     Biological homogeneity index    max     Datta & Datta, 2006
% 'BSI'     Biological stability index      max     Datta & Datta, 2006
% 'CH'      Calinski-Harabasz               max     Calinski & Harabasz, 1974
% 'CI'      C-index                         min     Hubert & Levin, 1976
% 'CON'     Connectivity index              min     Handl & Knowles, 2005
% 'DB'      Davies-Bouldin index            min     Davies & Bouldin, 1979
% 'DBMOD'   modified Davies-Bouldin index	min     Kim & Ramakrishna, 2005
% 'DN'      Dunn index                      max     Dunn, 1973
% 'DNG'     Dunn index using graphs         max     Pal & Biswas, 1997
% 'DNS'     modified Dunn index             max     Ilc, 2012
% 'FOM'     Figure of merit                 min     Yeung et al., 2001
% 'GAMMA'	Gamma index                     min     Baker & Hubert, 1975
% 'GDI_ij'  Generalized Dunn Indices (18)	max     Bezdek & Pal, 1998
% 'GPLUS'	G(+) index                      max     Rohlf, 1974
% 'HA' *    Hartigan index                  elbow	Hartigan, 1985
% 'HOM'     Homogeneity (average)           max     Sharan et al., 2003
% 'HOMMIN'  Homogeneity (minimum)           max     Sharan et al., 2003
% 'I'       I-index or PBM                  max     Maulik & Bandyopadhyay, 2002
% 'KL' *    Krzanowski-Lai index            max     Krzanowski & Lai, 1988
% 'SD' *    SD index                        min     Halkidi et al., 2000
% 'SDBW'	S_Dbw index                     min     Halkidi & Vazirgiannis, 2001
% 'SEP'     Separation (average)            min     Sharan et al., 2003
% 'SEPMAX'  Separation (maximum)            min     Sharan et al., 2003
% 'SF'      Score Function                  max     Saitta et al., 2008
% 'SIL'     Silhouette index (average)      max     Rousseeuw, 1987
% 'SSI'     Simple Structure Index          max     Dolnicar et al., 1999
% 'TAU'     Tau index                       min     Rohlf, 1974
% 'VAR'     Variance index                  min     Handl & Knowles, 2005b
% 'WR' *    Weighted inter-intra index      drop	Strehl, 2002
% 'WRP'*    Penalized WR index              drop	Strehl, 2002
% 'XB'      Xie-Beni index                  min     Xie & Beni, 1991
% 'XBMOD'   modified Xie-Beni index         min     Kim & Ramakrishna, 2005
% -------------------------------------------------------------------------
%
%   options.dtype: distance measure type (not all of the indices consider this):
%                      1 - Euclidean distance or
%                      2 - Pearson correlation or
%                      any other string identifier that is passed to pdist
%
%   options.NC  : vector that contains numbers of clusters, e.g., [2 3 4 5]
%                 length(NC) and number of labels sets in labels must
%                 match.
%
%	options.CON_L : number of nearest neighbors considered when calculating
%					CON index. If empty, default value of 5 is taken.
%
%   options.DB_p : parameter of the Minkowski distance
%                           p = 1: cityblock distance
%                           p = 2: Euclidean distance (default)
%   options.DB_q : parameter of dispersion measure of a cluster
%                           q = 1: average Euclidean measure (default)
%                           q = 2: standard deviation of the distance
%
%   options.GDI_interdist : (string/int) distance between two clusters
%                               1 - 'single' (original Dunn)
%                               2 - 'complete'
%                               3 - 'average'
%                               4 - 'centroid'
%                               5 - 'avg2cent'
%                               6 - 'hausdorff'
%                               - if empty, compute all of them
%   options.GDI_intradist : (string/int) cluster diameter:
%                               1 - 'complete' (original Dunn)
%                               2 - 'average'
%                               3 - 'avg2cent'
%                               - if empty, compute all of them
%
%   options.clMethod :  when using stability measures (APN, AD, ADM, BSI,
%                       FOM), additional labels with deletion have to be
%                       supplied; if not, they will be computed with method
%                       clMethod. See pplk_runClustererDel for details.
%
%   options.params :    parameters for clMethod
%
%   options.labelsDel : matrix [n X d X nNCs] with clustering labels when
%                       clustering data without one column. (i,j,k) contains
%                       clustering label for data sample i without j-th feature
%                       when clustering into NC(k) clusters.
%                       Required by stability measures (APN, AD, ADM, BSI, FOM).
%                       If not provided, labelsDel is computed with
%                       options.clMethod using options.params.
%
%   options.genenames : (required by BHI and BSI) cell of genenames, one
%                                                 for each row in data
%
%   options.annotations : (required by BHI and BSI)
%                         structure array with fields:
%                         .geneID         (gene identifier, e.g., Affymetrix probe ID)
%                         .functionID     (functional class identifier, e.g., GO ID)
%                         .aspect         (optional - when used with GO)
%                         .evidence       (optional - when used with GO)
%   options.GO_aspect : (required by BHI and BSI) which aspects to include:
%                               'BP' | 'CC' | 'MF' (or any subset of them);
%                               leave empty to compute all of them.
%   options.GO_evidence : (required by BHI and BSI) which evidence to ignore,
%                               i.e.: 'EXP', 'IDA', 'IEA',... or any subset;
%                               leave empty to compute all of them.
%   options.SD_alpha : (required by SD) weighting factor (should equals
%                      Dis(labels(:,cmax)), where cmax is the maximum
%                      number of input clusters). If empty set to number of
%                      clusters in labels.
%   options.SHOW :  logical, whether to plot indices (1) or not (0, default)
%
%   options.graph : already built graph on data for DNs and DNg
%   options.graph_type: type of a graph for DNs, DNg; default: 'gabriel'
%   options.uniqueInd : indeces of unique data points, required by DNs and
%                       DNg
%   options.graph_sqEucl: if weights of a graph are squared Euclidean
%                         distances, default 0
%   options.penalFun: penal function for increasing K (only for DNS) 
%                     ['none', 'reciproc', 'exponent', default: 'logistic']
%   options.penalFunP: strength of penals [0,1], 0: no penal, 1: full penal
%   options.penalFunLogMid : mid-point for logistic function, 
%                            default ceil(sqrt(N)/2)
%
% OUTPUTS
%	validInt	: a struct with fields that are named according to the used
%				  methods, i.e., DB if methods={'DB'}.
%	list		: numeric matrix [nMethods X nClusterings] with values of indices
%--------------------------------------------------------------------------
%   Requirements:
%       MATLAB Statistics Toolbox 
%       (Silhouette function, crosstab, nanmean, pdist, pdist2)
%   Acknowledgements:
%       Cluster Validity Analysis Platform (CVAP) (Version 3.4) 
%       Copyright (C) 2006-2007 by Kaijun Wang.
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- VERSION ----------------------------------------------------------
% Version: 1.2.1
% Last modified: 23-April-2014 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


callDir=chdir(pplk_homeDir());

[nrow, dim] = size(data);

if ~iscell(methods)
    methods = {methods};
end

%--------------------------------------------------------------------------
% OPTIONS check

% Defaults
dtype = 'euclidean';
genenames = [];
annotations = [];
GO_aspect = [];
GO_evidence = [];
clMethod = [];
clMethodParams = [];
CON_L=5;
DB_p = 2;
DB_q = 1;
GDI_interdist = 1:6;
GDI_intradist = 1:3;
GDI_mode = 'matrix'; % return num_inter X num_intra indices
I_p = 2;
SD_alpha = [];
NC = [];
SHOW = 0;
graph = [];
penalFun = 'logistic';
penalFunP = 1;
penalFunLogMid = ceil(sqrt(nrow)/2);

% Fill in with options supplied
if exist('options','var') && isstruct(options)
    fList = fieldnames(options);
    for fInd=1:length(fList)
        f= fList{fInd};
        switch lower(f)
            case 'dtype'
                dtype=options.(f);
            case 'genenames'
                genenames = options.(f);
            case 'annotations'
                annotations = options.(f);
            case 'go_aspect'
                GO_aspect = options.(f);
            case 'go_evidence'
                GO_evidence = options.(f);
            case 'clmethod'
                clMethod = options.(f);
            case 'clmethodparams'
                clMethodParams = options.(f);
            case 'db_p'
                DB_p = options.(f);
            case 'db_q'
                DB_q = options.(f);
            case 'gdi_interdist'
                if ~isempty(options.(f))
                    GDI_interdist = options.(f);
                end
            case 'gdi_intradist'
                if ~isempty(options.(f))
                    GDI_intradist = options.(f);
                end
            case 'i_p'
                I_p = options.(f);
            case 'con_l'
                CON_L = options.(f);
            case 'sd_alpha'
                SD_alpha = options.(f);
            case 'nc'
                NC = options.(f);
            case 'show'
                SHOW = options.(f);
            case 'graph'                
                graph = options.(f);
            case 'graph_type'                
                graph_type = options.(f);
            case 'uniqueind' 
                uniqueInd = options.(f);
            case 'graph_sqeucl'     
                graph_sqEucl = options.(f);
            case 'penalfun'     
                penalFun = options.(f);    
            case 'penalfunp'     
                penalFunP = options.(f);     
            case 'penalfunlogmid'     
                penalFunLogMid = options.(f);     
            otherwise
                warning('pplk:options',['Options field ',f,' is not relevant.']);
        end
    end
end

% Convert to uppercase
methods = upper(methods);

% Errors
if any(ismember(methods,{'BHI','BSI'}))
    if isempty(genenames)
        error('options.genenames is empty or does not exist! Gene names are required.');
    end
    if isempty(annotations)
        error('options.annotations is empty or does not exist! Gene annotations are required.');
    end
end

% Resolve GDI case - search for GDI_ij strings
GDImask = find(~cellfun(@isempty,strfind(methods,'GDI')));
GDIstr = methods(GDImask);
if ~isempty(GDImask)
    if any(strcmpi(GDIstr,'GDI')) && length(GDImask) > 1
        error('Do not mix GDI and GDI_ij names!');
    elseif ~any(strcmpi(GDIstr,'GDI'))
        if ~isequal(GDImask(1):GDImask(end), GDImask(1):GDImask(1)+length(GDImask)-1)
            error('GDI_ij names must be in sequence!');
        else
            GDI_interdist_s = [];
            GDI_intradist_s = [];
            for i = 1:length(GDIstr)
                S = GDIstr{i};
                if length(S) == 6
                    GDI_interdist_s = [GDI_interdist_s,str2num(S(5))];
                    GDI_intradist_s = [GDI_intradist_s,str2num(S(6))];
                end
            end
        end
        methods{GDImask(1)} = 'GDI';
        methods(GDImask(2:end)) = [];
        [GDI_interdist, ~,GDI_interLoc] = unique(GDI_interdist_s);
        [GDI_intradist, ~,GDI_intraLoc] = unique(GDI_intradist_s);
        GDI_mode = 'single'; % return only num_inter (or num_intra) indices
    end
end


%--------------------------------------------------------------------------

% Do we have multiple clusterings to validate?
% "Repair" labels if neccessary.
numLabels = size(labels,2);
if isempty(NC)
    if numLabels==1
        % ensure that labels are integers 1:K
        [clusterIDs, ~, labels] = unique(labels);
        NC = length(clusterIDs);        
    else
        NC = zeros(1,numLabels);
        for NCi = 1:numLabels
            % ensure that labels are integers 1:K
            [clusterIDs, ~, labels(:,NCi)] = unique(labels(:,NCi));
            NC(NCi) = length(clusterIDs);            
        end
    end
end

N = length(NC);

%--------------------------------------------------------------------------
% User defined functions
f=fopen(['..',filesep,'validation',filesep,'validInt.info'],'r');
usrfile=textscan(f,'%s %s %s %*[^\n]','delimiter','\t','commentStyle','%');
fclose(f);

usrMeth=usrfile{1};
usrAbbr=usrfile{2};
usrDesc=usrfile{3};

%--------------------------------------------------------------------------
% Identify meaningful groups of indices to speed up the computation

allMethods=[{
    'APN','AD','ADM',...
    'BHI','BSI',...
    'CH','CI','CON',...
    'DB','DB*', 'DN','DNG','DNS',...
    'FOM',...
    'GAMMA', 'GDI', 'GPLUS',...
    'HA','HOM','HOMMIN',...
    'I',...
    'KL',...
    'SD','SDBW','SEP','SEPMAX','SF','SIL','SSI',...
    'TAU',...
    'VAR',...
    'WR','WRP',...
    'XB','XB*'},...
    usrAbbr'];

if ~exist('methods','var') || isempty(methods) || all(strcmpi(methods,'all'))
    disp('Computing all CVI, incl. user-defined!');
    methods = allMethods;
end

nMethods=length(methods);

% determine, which group of methods will be launched
group_1 = {'CH','KL','HA'};
group_HOM = {'HOM','HOMMIN'};
group_SEP = {'SEP','SEPMAX'};
group_WR = {'WR','WRP'};
group_S = {'APN', 'AD', 'ADM', 'BSI', 'FOM'};
group_G = {'GAMMA', 'GPLUS', 'TAU'};
group_DNmod = {'DNG','DNS'};
group_SD = {'SD','SDBW'};
group_XB = {'XB','XB*','XBMOD'};
group_Dmat = [group_HOM, group_SEP, group_WR, group_G, ...
                {'APN','AD','ADM','FOM','CI','GDI','DN'}];
group_centers = {'DB','DB*','DBMOD','GDI','I','SSI','XB','XB*','XBMOD'};

groupMode1 = 0;
groupModeHOM = 0;
groupModeSEP = 0;
groupModeWR = 0;
groupModeS = 0;
groupModeG = 0;
groupModeDNmod = 0;
groupModeSD = 0;
groupModeXB = 0;
groupModeDmat = 0;
groupModeCenters = 0;

if any(ismember(group_1,methods))
    groupMode1 = 1;
end
if any(ismember(group_HOM,methods))
    groupModeHOM = 1;
end
if any(ismember(group_SEP,methods))
    groupModeSEP = 1;
end
if any(ismember(group_WR,methods))
    groupModeWR = 1;
end
if any(ismember(group_S,methods))
    groupModeS = 1;
end
if any(ismember(group_G,methods))
    groupModeG = 1;
end
if any(ismember(group_DNmod,methods))
    groupModeDNmod = 1;
end
if all(ismember(group_SD,methods))
    groupModeSD = 1;
end
if any(ismember(group_XB,methods))
    groupModeXB = 1;
end
if any(ismember(group_Dmat,methods))
    groupModeDmat = 1;
end
if any(ismember(group_centers,methods))
    groupModeCenters = 1;
end

% any of groups active? (ignore groupModeDmat and groupModeCenters)
groupModeAny = (groupMode1 | groupModeHOM | groupModeSEP | groupModeWR |...
                groupModeS | groupModeG | groupModeDNmod | groupModeSD | groupModeXB);

%--------------------------------------------------------------------------
% Calculating dis-similarity/distance matrix/vector of a data set
dname = dtype;
if (ischar(dtype) && (strcmpi(dtype,'pearson') || strcmpi(dtype,'correlation'))) || (isnumeric(dtype) && dtype==2)
    dname = 'correlation';
end

if groupModeDmat
    if isnumeric(data)
        if strcmp(dname,'correlation')
            % Pearson similarity [-1,1] is transformed to Pearson distance [0,1]             
            Dvector = pdist(data,'correlation')/2;
        else
            Dvector = pdist(data,dname);
        end
        Dmatrix = squareform(Dvector);
    else
        error('Input data should be numeric!');
    end
end

% Calculating cluster centers
if groupModeCenters
    sepMat = cumsum([0 NC]);
    tot = sum(sepMat);
    centers = zeros(tot,dim);
    U = false(tot,nrow);

    for i=1:N
        lbl = labels(:,i);
        E = eye(NC(i));
        Utmp= logical(E(:,lbl));
        U(sepMat(i)+1:sepMat(i+1),:) = Utmp;
        clsize = sum(Utmp,2);
        centers(sepMat(i)+1:sepMat(i+1),:) = bsxfun(@rdivide, Utmp*data, clsize);
    end    
end

%--------------------------------------------------------------------------
% Compute first the groups of indices
NCind = 1:N;

% skip if there is no groups
if groupModeAny
    
    % group_1 = {'CH','KL','HA'};
    if groupMode1
        CH = zeros(1,N);
        KL = zeros(1,N);
        HA = zeros(1,N);
    end
    
    % group_HOM _SEP _WR = {'HOM','HOMMIN','SEP','SEPMAX','WR','WRP'};
    if groupModeHOM
        HOM    = zeros(1,N);
        HOMMIN = zeros(1,N);
    end
    if groupModeSEP
        SEP    = zeros(1,N);
        SEPMAX = zeros(1,N);
    end
    if groupModeWR
        WR     = zeros(1,N);
        WRP    = zeros(1,N);
    end
    
    % group_S = {'APN','AD','ADM','BSI','FOM'};
    if groupModeS
        APN = zeros(1,N);
        AD = zeros(1,N);
        ADM = zeros(1,N);
        BSI = zeros(1,N);
        FOM = zeros(1,N);
    end
    
    % group_G = {'GAMMA','GPLUS','TAU'};
    if groupModeG
        GAMMA = zeros(1,N);
        GPLUS = zeros(1,N);
        TAU = zeros(1,N);
    end
    
    % group_DNmod = any{'DNG','DNS'};
    if groupModeDNmod
        graphOpt = [];
        
        % compute graph (only once) if it is not provided by options
        if isempty(graph)
            if ~exist('graph_type','var') || isempty(graph_type)
                graph_type = 'gabriel';
            end
            if ~exist('graph_sqEucl','var') || isempty(graph_sqEucl)
                graph_sqEucl = 0;
            end
            graphOpt.graph_sqEucl = graph_sqEucl;
            
            oldDir = chdir(pplk_homeDir(['validation',filesep,'indexDNmod']));            
            [graph,~,uniqueInd] = graph_create(data,[],graph_type,graphOpt);
            chdir(oldDir);
        % graph provided, check also provided parameters
        else            
            if ~exist('graph_sqEucl','var')
                error('Graph provided but graph_sqEucl missing!');
            end
            if ~exist('graph_type','var')
                error('Graph provided but graph_type missing!');
            end
            if ~exist('uniqueInd','var')
                error('Graph provided but uniqueInd missing!');
            end
            graphOpt.graph_sqEucl = graph_sqEucl;
        end
        % is graph directed?
        isDirected = strcmpi(graph_type, 'directedknn');
        
        % if the removal of duplicates occured, update data and labels
        if ~isempty(uniqueInd)
            dataGraphNew = data(uniqueInd,:);
            labelsGraphNew = labels(uniqueInd,:);
        else
            dataGraphNew = data;
            labelsGraphNew = labels;
        end
    end
    
    % group_SD = {'SD','SDBW'};
    if groupModeSD
        SD = zeros(1,N);
        SDBW = zeros(1,N);
    end
    
    % group_XB = {'XB','XB*'};
    if groupModeXB
        XB = zeros(1,N);
        XBMOD = zeros(1,N);
    end
    
    % loop over k
    for k = NCind
        
        if groupMode1
            % Davies-Bouldin, Calinski-Harabasz, Dunn, Krzanowski-Lai, and Hartigan
            % index - be aware that KL and HA work only for sequence of clusterings
            % with different number of clusters (e.g., from 2 to 10).
            
            % Dunn is probably returning wrong results (using external function)!
            % Davies-Bouldin is implemented in external file and gives identical
            % results as this one.
            [DBdummy,CH(k),DNdummy,KL(k),HA(k),ST] = valid_internal_deviation(data,labels(:,k),dtype);
        end
        
        if groupModeHOM || groupModeSEP || groupModeWR
            oldDir = chdir(['..' filesep 'validation' filesep 'indexHOM_SEP_WR']);
            
            if groupModeWR
                [HOM_tmp,HOMMIN_tmp,SEP_tmp,SEPMAX_tmp,WR_tmp,WRP_tmp] = ...
                        indexHOM_SEP_WR(Dmatrix,labels(:,k),'dist',dtype);
            else
                if groupModeSEP
                    [HOM_tmp,HOMMIN_tmp,SEP_tmp,SEPMAX_tmp] = ...
                        indexHOM_SEP_WR(Dmatrix,labels(:,k),'dist',dtype);
                else
                    [HOM_tmp,HOMMIN_tmp] = ...
                        indexHOM_SEP_WR(Dmatrix,labels(:,k),'dist',dtype);
                end
            end
            
            if groupModeHOM
                HOM(k) = HOM_tmp;
                HOMMIN(k) = HOMMIN_tmp;
            end
            if groupModeSEP
                SEP(k) = SEP_tmp;
                SEPMAX(k) = SEPMAX_tmp;
            end
            if groupModeWR
                WR(k) = WR_tmp;
                WRP(k) = WRP_tmp;
            end
            chdir(oldDir);            
        end
        %if groupMode2
            % weighted inter/intra ratio
            % S = ind2cluster(labels(:,k));
            
            % the last parameter is set to 0 in order to compute all indices. To
            % compute only CI, change it to 1
            % [HOM(k), SEP(k), dummyCI, WR(k)] = valid_internal_intra(Dmatrix, S, dtype, 0);
        %end
        
        if groupModeS
            
            % check if labels from clustering with deletion are available. If
            % not, create them using options.clMethod fed with options.clMethodParams
            if exist('options','var') && isstruct(options) && isfield(options,'labelsDel')
                labelsDel = options.labelsDel(:,:,k);
            else
                if ~isempty(clMethod)
                    labelsDel = pplk_runClustererDel(clMethod,data,NC(k),1,clMethodParams);
                    
                else
                    error('Running stability method without labelsDel provided and clustering method unspecified. Provide with options.clMethod.');
                end
                
            end
            % BSI
            if any(strcmpi('BSI',methods))
                oldDir = chdir(['..' filesep 'validation' filesep 'biological']);
                BSI(k) = indexBSI(genenames, labels(:,k), labelsDel, annotations, GO_aspect, GO_evidence);
                chdir(oldDir);
            end
            
            % stability measures
            if any(ismember(methods,{'APN','AD','ADM','FOM'}))
                oldDir = chdir(['..' filesep 'validation' filesep 'stability']);
                
                [APN(k), AD(k), ADM(k), FOM(k)] = ...
                    stability(data, Dmatrix, labels(:,k), labelsDel, dtype);
                
                chdir(oldDir);
            end
        end
        
        if groupModeG
            oldDir = chdir(['..' filesep 'validation' filesep 'indexGammaGplusTau']);
            [G, Gmod, GPLUS(k), TAU(k)] = indexGammaGplusTau(Dvector,labels(:,k),1);
            GAMMA(k) = G; % = Gmod % Gmod is a normalized version
            chdir(oldDir);
        end
        
                
        % If both SD and SDBW are selected
        if groupModeSD
            oldDir = chdir(['..' filesep 'validation' filesep 'indexSD_SDbw']);
            
            if isempty(SD_alpha)
                [~,cmaxi] = max(NC); % index of labels with max. num. of clust.
                SD_alpha = SD_Dis(data,labels(:,cmaxi));
            end
            
            [SD(k), SDBW(k)] = indexSD_SDbw(data,labels(:,k),SD_alpha);
            chdir(oldDir);
        end
        
        if groupModeXB
            oldDir = chdir(['..' filesep 'validation' filesep 'indexXB']);
            lblStruct.U = U(sepMat(k)+1:sepMat(k+1),:);
            lblStruct.center = centers(sepMat(k)+1:sepMat(k+1),:); 
            [XB(k), XBMOD(k)] = indexXB(data, lblStruct);
            chdir(oldDir);
        end
    end
    
    
    if groupMode1
        % If there are multiple clusterings to evaluate (to find number of
        % clusters), some adjustments have to be performed
        if N > 1
            S = trace(ST);
            KL = [S KL];
            HA = [S HA];
            R = abs(KL(1:N)-KL(2:N+1));
            S = [R(2: end) R(end)];
            KL = R./S;
            KL(N) = KL(N-1);
            R = HA(1:N)./HA(2:N+1);
            HA = (R-1).*(nrow-[NC(1)-1 NC(1:N-1)]-1);
        end
    end
end

%--------------------------------------------------------------------------
% Store the results (of already computed group indices) and compute all the
% others indices

% if we have to compute GDI indices, increase number of methods accordingly
if (any(strcmpi('GDI',methods)))
    switch GDI_mode
        case 'matrix'
            nGDI = length(GDI_interdist)*length(GDI_intradist);
        case 'single'
            nGDI = length(GDI_interLoc);
    end
else
    nGDI = 1;
end
% reserve space - if GDI, add number of variants (-1 because 'GDI' is
% already in nMethods)
list = zeros(nMethods+nGDI-1,N);
ind = 1;

for mInd=1:nMethods
    currMethod=upper(methods{mInd});
    
    switch(currMethod)
        
        case 'AD'
            validInt.AD=AD;
            list(ind,:)=AD; ind = ind+1;
            
        case 'ADM'
            validInt.ADM=ADM;
            list(ind,:)=ADM; ind = ind+1;
            
        case 'APN'
            validInt.APN=APN;
            list(ind,:)=APN; ind = ind+1;
            
        case 'BHI'
            oldDir = chdir(['..' filesep 'validation' filesep 'biological']);
            BHI = zeros(1,N);
            for k = NCind
                BHI(k) = indexBHI(genenames, labels(:,k), annotations, GO_aspect, GO_evidence);
            end
            validInt.BHI=BHI;
            list(ind,:)=BHI; ind = ind+1;
            chdir(oldDir);
            
        case 'BSI'
            validInt.BSI=BSI;
            list(ind,:)=BSI; ind = ind+1;
            
        case 'CH'
            validInt.CH=CH;
            list(ind,:)=CH; ind = ind+1;
            
        case 'CI'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexCI']);
            CI = zeros(1,N);
            for k = NCind
                CI(k) = indexCI(Dvector, labels(:,k),1);
            end
            validInt.CI=CI;
            list(ind,:)=CI; ind = ind+1;
            chdir(oldDir);
            
        case 'CON'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexCON']);
            CON = zeros(1,N);
            for k = NCind
                CON(k) = indexCON(data,labels(:,k),CON_L);
            end
            validInt.CON=CON;
            list(ind,:)=CON; ind = ind+1;
            chdir(oldDir);
            
        case 'DB'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexDB']);
            DB = zeros(1,N);
            for k = NCind
                DB(k) = indexDB(data, labels(:,k),0,DB_p,DB_q,centers(sepMat(k)+1:sepMat(k+1),:));
            end
            validInt.DB = DB;
            list(ind,:)=DB; ind = ind+1;
            chdir(oldDir);
            
        case {'DBMOD','DB*'}
            oldDir = chdir(['..' filesep 'validation' filesep 'indexDB']);
            DBMOD = zeros(1,N);
            for k = NCind
                DBMOD(k) = indexDB(data, labels(:,k),1,DB_p,DB_q,centers(sepMat(k)+1:sepMat(k+1),:));
            end
            validInt.DBMOD = DBMOD;
            list(ind,:)=DBMOD; ind = ind+1;
            chdir(oldDir);
            
        case 'DN'
            % using external method - works fine
            oldDir = chdir(['..' filesep 'validation' filesep 'indexDN']);
            DN = zeros(1,N);
            for k = NCind
                DN(k) = dunns(NC(k),Dmatrix, labels(:,k));
            end
            validInt.DN=DN;
            list(ind,:)=DN; ind = ind+1;
            chdir(oldDir);
            
        case {'DNG'}
            oldDir = chdir(['..' filesep 'validation' filesep 'indexDNmod']);
            DNG = zeros(1,N);
            for k = NCind
                %repair labels - removal of duplicates in
                % data may cause inconsistency in labels (if singleton cluster
                % was removed, one integer is missing in labels)
                [~,~,labelsNewRepaired] = unique(labelsGraphNew(:,k));
                DNG(k) = indexDNg_graph(graph,labelsNewRepaired,dataGraphNew);
            end
            validInt.DNG = DNG;
            list(ind,:)=DNG; ind = ind+1;
            chdir(oldDir);
            
        case {'DNS'}            
            oldDir = chdir(['..' filesep 'validation' filesep 'indexDNmod']);
            DNS = zeros(1,N);
            for k = NCind
                %repair labels - removal of duplicates in
                % data may cause inconsistency in labels (if singleton cluster
                % was removed, one integer is missing in labels)
                [~,~,labelsNewRepaired] = unique(labelsGraphNew(:,k));
                DNS(k) = indexDNs_graph(graph,dataGraphNew,labelsNewRepaired,isDirected,graph_type,graphOpt,penalFun,penalFunP,penalFunLogMid);
            end
            validInt.DNS = DNS;
            list(ind,:) = DNS; ind = ind+1;
            chdir(oldDir);
            
            
        case 'FOM'
            validInt.FOM=FOM;
            list(ind,:)=FOM; ind = ind+1;
            
        case 'GAMMA'
            validInt.GAMMA=GAMMA;
            list(ind,:)=GAMMA; ind = ind+1;
            
        case 'GDI'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexGDI']);
            GDI = zeros(nGDI,N);
            for k = NCind
                GDImat = indexGDI(data,labels(:,k),GDI_interdist,GDI_intradist,Dmatrix,centers(sepMat(k)+1:sepMat(k+1),:),dtype,dtype);
                
                switch GDI_mode
                    case 'matrix'
                        GDImat=GDImat';
                        GDI(:,k)=GDImat(:);
                    case 'single'
                        GDI(:,k) = GDImat(sub2ind(size(GDImat),GDI_interLoc,GDI_intraLoc));
                end
            end
            
            for ii = 1:length(GDIstr)
                validInt.(GDIstr{ii})=GDI(ii);
            end
            
            list(ind:ind+nGDI-1,:) = GDI; ind = ind+nGDI;
            chdir(oldDir);
            
        case 'GPLUS'
            validInt.GPLUS=GPLUS;
            list(ind,:)=GPLUS; ind = ind+1;
            
        case 'HA'
            validInt.HA=HA;
            list(ind,:)=HA; ind = ind+1;
            
        case 'HOM'
            validInt.HOM=HOM;
            list(ind,:)=HOM; ind = ind+1;
        case 'HOMMIN'
            validInt.HOMMIN=HOMMIN;
            list(ind,:)=HOMMIN; ind = ind+1;
            
        case 'I'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexI']);
            I = zeros(1,N);
            for k = NCind
                I(k) = indexI(data, labels(:,k),I_p,centers(sepMat(k)+1:sepMat(k+1),:));                
            end
            validInt.I=I;
            list(ind,:)=I; ind = ind+1;
            chdir(oldDir);
            
        case 'KL'
            validInt.KL=KL;
            list(ind,:)=KL; ind = ind+1;
            
        case 'SEP'
            validInt.SEP=SEP;
            list(ind,:)=SEP; ind = ind+1;
        case 'SEPMAX'
            validInt.SEPMAX=SEPMAX;
            list(ind,:)=SEPMAX; ind = ind+1;
            
        case 'SD'
            if groupModeSD
                validInt.SD=SD;
                list(ind,:)=SD; ind = ind+1;
            else
                oldDir = chdir(['..' filesep 'validation' filesep 'indexSD_SDbw']);
                if isempty(SD_alpha)
                    [~,cmaxi] = max(NC); % index of labels with max. num. of clust.
                    SD_alpha = SD_Dis(data,labels(:,cmaxi));
                end
                SD = zeros(1,N);
                for k = NCind
                    SD(k) = indexSD(data, labels(:,k),SD_alpha);
                end
                validInt.SD=SD;
                list(ind,:)=SD; ind = ind+1;
                chdir(oldDir);
            end
            
        case 'SDBW'
            if groupModeSD
                validInt.SDBW=SDBW;
                list(ind,:)=SDBW; ind = ind+1;
            else
                oldDir = chdir(['..' filesep 'validation' filesep 'indexSD_SDbw']);
                SDBW = zeros(1,N);
                for k = NCind
                    SDBW(k) = indexSDbw(data, labels(:,k));
                end
                validInt.SDBW=SDBW;
                list(ind,:)=SDBW; ind = ind+1;
                chdir(oldDir);
            end
            
        case 'SF'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexSF']);
            SF = zeros(1,N);
            for k = NCind
                SF(k) = indexSF(data, labels(:,k));
            end
            validInt.SF=SF;
            list(ind,:)=SF; ind = ind+1;
            chdir(oldDir);
            
        case 'SIL'
            % Silhouette index - MATLAB implementation (Statistics toolbox)
            SIL = zeros(1,N);
            for k = NCind
                s = silhouette(data, labels(:,k), dname);
                SIL(k) = mean(s);
            end
            validInt.SIL=SIL;
            list(ind,:)=SIL; ind = ind+1;
            
        case 'SSI'    
            oldDir = chdir(['..' filesep 'validation' filesep 'indexSSI']);
            SSI = zeros(1,N);
            for k = NCind
                SSI(k) = indexSSI(centers(sepMat(k)+1:sepMat(k+1),:), labels(:,k), 1);            
            end
            validInt.SSI = SSI;
            list(ind,:)=SSI; ind = ind+1;
            chdir(oldDir);
            
        case 'TAU'
            validInt.TAU=TAU;
            list(ind,:)=TAU; ind = ind+1;
            
        case 'VAR'
            oldDir = chdir(['..' filesep 'validation' filesep 'indexVAR']);
            VAR = zeros(1,N);
            for k = NCind
                VAR(k) = indexVAR(data, labels(:,k));
            end
            validInt.VAR = VAR;
            list(ind,:)=VAR; ind = ind+1;
            chdir(oldDir);
            
        case 'WR'
            validInt.WR=WR;
            list(ind,:)=WR; ind = ind+1;
        case 'WRP'
            validInt.WRP=WRP;
            list(ind,:)=WRP; ind = ind+1;
            
        case 'XB'
            validInt.XB = XB;
            list(ind,:) = XB; ind = ind+1;
            
        case {'XBMOD','XB*'}
            validInt.XBMOD = XBMOD;
            list(ind,:) = XBMOD;  ind = ind+1;
            
        otherwise
            % Check for user-defined imported functions, relating file
            % validInt.info .
            
            % search is made in the list of abbreviations, stored in
            % usrAbbr string array.
            idx=strcmp(currMethod,upper(usrAbbr));
            
            if any(idx)
                disp(['USER-DEFINED :: Executing function: ',usrAbbr{idx},'=',usrMeth{idx},'(...) -- ',usrDesc{idx}])
                oldPath=chdir(['..',filesep,'validation']);
                USER=zeros(1,N);
                for k=NCind
                    USER(k)=feval(str2func(usrMeth{idx}),data,labels(:,k));
                end
                validInt.(usrAbbr{idx})=USER;
                list(ind,:) = USER; ind = ind+1;
                chdir(oldPath);
            else
                disp(['Non-existing method name: ', methods{mInd}, '! Ignoring.']);
                validInt.(methods{mInd}) = NaN;
                list(ind,:) = NaN; ind = ind+1;
            end
    end
end
%--------------------------------------------------------------------------
% Plot results if reqired
if SHOW
    f=fieldnames(validInt);
    for i=f(:)'
        indVal=validInt.(i{1});
        figure();
        plot(indVal,'-*');
        title(i);
    end
end

%--------------------------------------------------------------------------
% Change back to the calling folder
chdir(callDir);

end

%=========================================================================


function [clusters, newlabels] = ind2cluster(labels)
% Transform vector of integers (cluster labels) to cell array.
C = unique(labels);
newlabels = labels;
k = length(C);
clusters = cell(1,k);
for i = 1:k
    ind = find(labels==C(i));
    clusters{i} = ind;
    newlabels(ind) = i;
end

end


function [DB,CH,Dunn,KL,Han,st] = valid_internal_deviation(data,labels,dtype)
% cluster validity indices based on deviation & sum of squares

[nrow,nc] = size(data);
labels = double(labels);
k = max(labels);
if strcmp(dtype,'euclidean') || (isnumeric(dtype) && dtype == 1)
    [st,sw,sb,cintra,cinter] = valid_sumsqures(data,labels,k);
else
    [st,sw,sb,cintra,cinter] = valid_sumpearson(data,labels,k);
end
ssw = trace(sw);
ssb = trace(sb);

if k > 1
    % Davies-Bouldin & Dunn based on centroid diameter & linkage distance
    %[DB, Dunn] = valid_DbDunn(cintra, cinter, k);
    DB = NaN;
    Dunn = NaN;
    CH = ssb/(k-1);
else
    CH =ssb;
    DB = NaN;
    Dunn = NaN;
end

CH = (nrow-k)*CH/ssw;    % Calinski-Harabasz
Han = ssw;               % component of Hartigan
KL = (k^(2/nc))*ssw;     % component of Krzanowski-Lai

end

function [T, W, B, Sintra, Sinter] = valid_sumsqures(data,labels,k)
%   data:    a matrix with each column representing a variable.
%   labels: a vector indicating class labels
%   W         within-group sum of squares and cross-products
%   B           between-group sum of squares and cross-products
%   T           total sum of squares and cross-products
% Sintra & Sinter: centroid diameter & linkage distance

if (size(labels, 1) == 1)
    labels = labels';
end

[ncase,m] = size(data);

% computing the Total sum of squares matrix
Dm = mean(data);
Dm = data - Dm(ones(ncase,1),:);
T = Dm'*Dm;

% computing within sum of squares matrix
W = zeros(size(T));
Dm = zeros(k,m);
Sintra = zeros(1,k);
for i = 1:k
    if k > 1
        Cindex = find(labels == i);
    else
        Cindex = 1:ncase;
    end
    nk = length(Cindex);
    if nk > 1
        dataC = data(Cindex,:);
        m = mean(dataC);
        Dm(i,:) = m;
        dataC = dataC - repmat(m,nk,1);  %m(ones(nk,1),:)
        W = W + dataC'*dataC;
        dataC = sum(dataC.^2,2);
        Sintra(i) = mean(sqrt(dataC));  % distances to cluster center
    end
end

B = T - W;

% distances between cluster centers
Sinter = zeros(k,k);
if k > 1
    for i = 1:k
        for j = i+1:k
            m = abs(Dm(i,:) - Dm(j,:));
            Sinter(i,j) = sqrt(sum(m.^2));
            Sinter(j,i) = Sinter(i,j);
        end
    end
end

end

function [T, W, B, Sintra, Sinter] = valid_sumpearson(data,labels,k)
% within-, between-cluster and total sum of squares
% Sintra/Sinter: centroid diameter/linkage distance based on Pearson correlation

C = mean(data);
R = similarity_pearsonC(data', C');
T = R*R';

W = 0;
Sintra = zeros(1,k);
Sinter = zeros(k,k);
for i = 1:k
    Ui = find(labels == i);
    ni = length(Ui);
    if ni > 1
        datai = data(Ui,:);
        C = mean(datai);
        R = similarity_pearsonC(datai', C'); % distances to cluster center
        Sintra(i) = mean(R);
        W = W + R*R';
    end
    % distances between cluster centers
    for j = i+1:k
        Ui = find(labels == j);
        ni = length(Ui);
        if ni > 0
            datai = data(Ui,:);
            if ni == 1
                Cj = datai;
            else
                Cj = mean(datai);
            end
            Sinter(i,j) = similarity_pearsonC(Cj', C');
        end
    end
end

B = T - W;

end

function R = similarity_pearsonC(data, C)
% pearson coefficients between every column and the center
% input matrix: data --- nrow rows * ncol columns
% output matrix:   R --- ncol columns

[nrow,ncol] = size(data);
dm = mean(data);
data = data-repmat(dm,nrow,1);
C = C-mean(C);
R = ones(1,ncol);

X = sqrt(C'*C);
for j = 1:ncol
    y = data(:,j);
    xy = C'*y;
    Y = sqrt(y'*y);
    S = X*Y;
    % if S == 0      S = NaN;   end
    R(j) = xy/S;
end

% Pearson similarity [-1,1] is normalized to Pearson distance [0,1]
R = 1-(1+R)*0.5;
end