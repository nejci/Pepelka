function [labelsEns, moreInfo] = pplk_genEns(data, methods, params)
% [labelsEns, moreInfo] = pplk_genEns(data, methods, params)
% Function pplk_genEns creates ensemble, according to methods cell. It
% performs data clustering with selected algoritms and packs the results
% into a N-by-ensembleSize labelsEns matrix.
%
% INPUTS		
%   data	
%       A n-by-d matrix of input data.
%
%   methods		
%       A M-by-4 cell, containing details about ensemble structure:
%       {methodName1, repetitions1, k1, mode1; ...
%       methodName2, repetitions2, k2, mode2; ...
%       ...
%       methodNameM, repetitionsM, kM, modeM}
%
%   params			
%       Parameters structure; can be [] or non-existent for defaults;
%       params.subsampling controls data subsampling before clustering
%       subsmpl = params.subsampling is a cell with 1 or 2 elements:
%           subsmpl{1} - type:
%               'none'
%                   No subsampling.
%               'rows'
%                   Sample by rows.
%               'cols'
%                   Sample by columns
%
%           subsmpl{2} - size of sample:
%               percentage
%                   Number on interval (0,1)); 0.1 means 10%.
%               percentage interval [a,b]
%                   Choose randomly from interval. Number on interval 
%                   [1, no.elements].
%               string 'rand'
%                   Random number on interval [1, no.elements].
%
%
% OUTPUTS		
%   labelsEns
%       Clusterings - N-by-1 column vectors, obtained from single
%       clusterer, concatenated into a N-by ensembleSize matrix.
%
%   moreInfo
%       Other info from clusterer (time, etc.).
%
%
% DETAILS on inputs 'methods' and 'params'
%
% List of available methods with associated parameters:
%   Legend:
%   methodName (params.methodName_[val]) Description
%   *  can not-exists
%   ** can be []
%
%   - AL (distance*) Average-linkage.
%   - CS (sigma**, Kin, Nin ) Cauchy-Schwarz divergence.
%   - EM (/) Expectation Maximization.
%   - FCM (/) Fuzzy C-means.
%   - GSOM (G,dG,alfa,maxIter, pOut, msize**, shape, distance, showSOM, 
%     showGrav, advanced) Gravitational SOM.        
%   - HCL (clustMethod, distance) Hierarchical Clustering (general,
%     including AL,SL,WL)
%   - SOMKM (nRuns, msize**, shape) K-means on the SOM.
%   - KCC (/) K-centers.
%   - KM (maxIter, nRuns, distance) K-means.
%   - KVV (sigma*) Kannan, Vempala and Vetta spectral   
%   - NC (/) Normalized cuts.      
%   - NJW (sigma*) Ng, Jordan and Weiss spectral.
%   - RANDOM (/) Random partition.
%   - SL (distance*) Single-linkage.
%   - WL (distance*) Ward-linkage 
%
%
% FORMAT of methods cell row:
%   methodName	
%       Select one from the list above.
%
%   repetitions 	
%       Number of repetitions of each method. Partition from each
%       repetition is regarded as 'ensemble member'. If mode is 'fixed' and
%       k is a vector [kmin,kmax], there will be
%       (repetitions*(kmax-kmin+1)) runs.
%
%   k
%       Desired number of clusters if mode='fixed' or upper bound for
%       number of clusters if mode='random'. Also, k can be a
%       vector of two elements [kmin,kmax]. If k is an empty matrix [], it
%       is automatically determined (only if clusterer
%       supports this).
%
%   mode 
%       'fixed'
%           k is target number (interval) of clusters.
%       'rand'  
%           Number of clusters is randomly chosen on [2,k] if k is a 
%           scalar, or on [kmin,kmax] if k is a vector.
%
%
% EXAMPLES using methods formats
%
%   1. Homogeneous ensemble:
%       a) K-means (100 runs) with fixed number of clusters (k = 3)
%           methods = {'KM',100,3,'fixed'};
%
%       b) K-means (100 runs) with fixed number of clusters (k = sqrt(N))
%           methods = {'KM',100,[],'fixed'};
%
%       c) K-means with fixed number of clusters: go from 10 to 20
%           methods = {'KM',1,[10,20],'fixed'};
%
%       d) K-means with fixed number of clusters: go from 10 to 20 and on
%          each step generate 2 partitions (repetitions)
%           methods = {'KM',2,[10,20],'fixed'};
%
%       e) K-means (100 runs) with random number of clusters (from 2 to 
%       sqrt(N))
%           methods = {'KM',100,[],'rand'};
%
%       f) K-means (100 runs) with random number of clusters (from 2 to 10)
%           methods = {'KM',100,10,'rand'};
%
%       g) K-means (100 runs) with random number of clusters (from 10 to 
%       20)
%           methods = {'KM',100,[10,20],'rand'};
%
%	2. Heterogeneous ensemble - based on various methods with fixed number
%	of clusters (k=10)
%       methods =   {'KM',  1,10,'fixed'; ...
%                   'GSOM',1,10,'fixed';...
%                   'NC',  1,10,'fixed'};
%
%	3. Mixed ensemble - based on various methods with various number of
%	runs and modes of choosing number of clusters - the most general case.
%       methods =   {'KM',10,10,'rand'; ...
%                   'GSOM',1,[],'fixed';...
%                   'NC',10,[5,10],'rand';...
%                   'NC',5,10,'fixed';...
%                   'AL',1,[5,10],'fixed'};
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka
% See also runClusterer, setParamsDefault, setParamsForData


callDir=chdir(pplk_homeDir());


% PARAMETERS - if params doesn't exist or is empty, load default values
%-------------------------------------------------------------------------

if ~exist('params','var') || isempty(params)
    params=pplk_setParamsDefault();
end


% RUN CLUSTERERS - generating cluster ensemble members, using different modes
%-------------------------------------------------------------------------
%		- 'fixed' : fixed number of clusters k or [k(1),k(2)]; default is k=sqrt(N)
%		- 'rand'  : pick a random number on interval [2,k] or [k(1),k(2)]
%
% TODO:
%	user methods integration (thru mode='userDefinedMode')


methods_name = upper(methods(:,1));
methods_R = cell2mat(methods(:,2));
methods_k = methods(:,3);
methods_mode = lower(methods(:,4));

% determine data matrix size
N=size(data,1);
dataOrig = data;

% Check whether there is subsampling required.
doSubsmpl = 0;
if isfield(params,'subsampling')
    % subsmpl = params.subsampling is a cell with 1 or 2 elements:
    % subsmpl{1} - type:
    %                    'none': no subsampling
    %                    'rows': sample by rows
    %                    'cols': sample by columns
    % subsmpl{2} - size of sample:
    %                    percentage; number on interval (0,1)); 0.1 means 10%
    %                    number on interval [1, no.elements]
    %                    string 'rand'; random number on interval [1, no.elements]
    
    subsmpl = params.subsampling;
    subsmpl_type = 'none';
    subsmpl_size = 'rand';
    
    if length(subsmpl) == 1
        subsmpl_type = subsmpl{1};
    elseif length(subsmpl) == 2
        subsmpl_type = subsmpl{1};
        subsmpl_size = subsmpl{2};
    else
        error('Wrong subsampling method!');
    end
    
    if ~strcmpi(subsmpl_type,'none')
        doSubsmpl = 1;
        
        if strcmpi(subsmpl_type,'rows')
            sub_dim = 1;
        elseif strcmpi(subsmpl_type,'cols')
            sub_dim = 2;
        end
        
        sub_size = get_sub_size(size(data,sub_dim),subsmpl_size);       
        
        if sub_dim == 1
            N = sub_size;
        end
    end
end

% determine for which row of methods k is a vector and mode is 'fixed'
k_length = cellfun(@length,methods_k);
fixed_mode_mask = strcmpi('fixed',methods_mode);
k_interval_mask = k_length == 2;
mask = k_interval_mask & fixed_mode_mask;
% Calculate ensemble size; consider interval span of k and multiply it by
% number of repetitions.
k_masked = cell2mat(methods_k(mask));
R_masked = methods_R(mask);
R_size = methods_R;
if ~isempty(k_masked)
    R_size(mask) = (k_masked(:,2)-k_masked(:,1)+1).*R_masked;
end

M = length(methods_name);
ensembleSize = sum(R_size);

labelsEns = zeros(N,ensembleSize);
moreInfo = cell(1,ensembleSize);  % additional info from clusterers (elapsed time, etc.)

%index of current position in matrix labels - tells which ensemble member
%is being considered
memberInd=1;

% loop goes through all methods in the list
for ind=1:M
    
    method = methods_name{ind};
    R = methods_R(ind);
    k = methods_k{ind};
    mode = methods_mode{ind};
    
    % default value for number of clusters in ensemble
%     if isempty(k)
%         k = ceil(sqrt(N)); %rule of thumb
%     end

    % special case, when method is SPECLS and there is no features
    % subsampling: precompute eigen vectors for speed
    data = dataOrig;

    if strcmpi(method,'SPECLS2') 
        if ~doSubsmpl
            Knn = 7;
            if isfield(params,'SPECLS_Knn')
                Knn = params.SPECLS_Knn;
            end
            SPECLS_dir = pplk_homeDir(['methods',filesep,'SPECLS']);
            oldDir = chdir(SPECLS_dir);
            data = SPECLS_eigv(data,Knn);
            chdir(oldDir);            
        else
            % feature subsampling not suported by SPECLS2
            method = 'SPECLS';
        end
    end
    if strcmpi(method,'SPECLS3') 
        if ~doSubsmpl
            % load precomputed SVD data
            Knn = 7;
            if isfield(params,'SPECLS_Knn')
                Knn = params.SPECLS_Knn;
            end
            dataName = params.dataName;           
            
            SPECLS_dir = pplk_homeDir(['methods',filesep,'SPECLS']);
            oldDir = chdir(SPECLS_dir);
            
            preV = load('precomputedV.mat','Knn','datasets');
            if any(preV.Knn == Knn) && any(strcmp(dataName,preV.datasets))
                V_file = ['V_',dataName,'_knn',num2str(Knn)];
                preV = load('precomputedV.mat',V_file);
                data = preV.(V_file);
            else
                data = SPECLS_eigv(data,Knn);
            end
            chdir(oldDir);            
        else
            % feature subsampling not suported by SPECLS3 - compute SVD for
            % each replications
            method = 'SPECLS3SUB';
        end
    end
    
    switch mode
        %run selected method R-times with fixed number of clusters k
        case {'fix','fixed'}
            if k_interval_mask(ind)
                kmin = k(1);
                kmax = k(2);
                
                for ki = kmin:kmax
                    for r=1:R
                        if doSubsmpl
                            sub_ind = sort(randsample(size(dataOrig,sub_dim),sub_size));
                            if sub_dim == 1
                                data = dataOrig(sub_ind,:);
                            else
                                data = dataOrig(:,sub_ind);
                                % randomize again
                                sub_size = get_sub_size(size(data,2),subsmpl_size); 
                            end
                        end
                        
                        try
                            [labels, mI] = pplk_runClusterer(method,data,ki,1,params);
                        catch err
                            labels = nan(N,1);
                            mI{1}.ERROR = err.message;
                        end
                        
                        mI{1}.methodsInfo={method,R,ki,mode};
                        if doSubsmpl
                            mI{1}.subsample = subsmpl;
                            mI{1}.subsampleInd = sub_ind;
                        end
                        
                        moreInfo{memberInd}= mI{1};
                        labelsEns(:,memberInd)= labels;
                        
                        memberInd=memberInd+1;
                    end
                end
            else
                for r=1:R
                    if doSubsmpl
                        sub_ind = sort(randsample(size(dataOrig,sub_dim),sub_size));
                        if sub_dim == 1
                            data = dataOrig(sub_ind,:);
                        else
                            data = dataOrig(:,sub_ind);
                            % randomize again
                            sub_size = get_sub_size(size(data,2),subsmpl_size);
                        end
                    end
                    
                    try
                        [labels, mI] = pplk_runClusterer(method,data,k,1,params);
                    catch err
                        labels = nan(N,1);
                        mI{1}.ERROR = err.message;
                    end
                    
                    mI{1}.methodsInfo={method,R,k,mode};
                    if doSubsmpl
                        mI{1}.subsample = subsmpl;
                        mI{1}.subsampleInd = sub_ind;
                    end
                    
                    moreInfo{memberInd}= mI{1};
                    labelsEns(:,memberInd)= labels;
                    
                    memberInd=memberInd+1;
                end
            end
            
            % run selected method R-times with random number of clusters
        case {'rand','random'}
            if k_interval_mask(ind)
                kmin = k(1);
                kmax = k(2);
            else
                kmin = 2;
                kmax = k;
            end
            
            % create R randomized numbers of clusters from [kmin,kmax]
            k_rand = randi([kmin,kmax],1,R);
            
            for r=1:R
                if doSubsmpl
                    sub_ind = sort(randsample(size(dataOrig,sub_dim),sub_size));
                    if sub_dim == 1
                        data = dataOrig(sub_ind,:);
                    else
                        data = dataOrig(:,sub_ind);
                        % randomize again
                        sub_size = get_sub_size(size(data,2),subsmpl_size);
                    end
                end
                
                
                try
                    [labels, mI] = pplk_runClusterer(method,data,k_rand(r),1,params);
                catch err
                    labels = nan(N,1);                    
                    mI{1}.ERROR = err.message;
                end
                
                mI{1}.methodsInfo={method,R,k,mode};
                mI{1}.methodsInfoKtrue=k_rand(r);
                if doSubsmpl
                    mI{1}.subsample = subsmpl;
                    mI{1}.subsampleInd = sub_ind;
                end
                
                moreInfo{memberInd} = mI{1};
                labelsEns(:,memberInd)= labels;
                
                memberInd=memberInd+1;
            end
            
        otherwise
            error(['Mode ',mode,' does not exist! Use ''fixed'' or ''rand''.']);
    end
end

chdir(callDir);

end

function sub_size = get_sub_size(data_dim_size,subsmpl_size)
% determine subset size of data
% data_dim_size = size(data,selected_dim);

if ischar(subsmpl_size)
    if strcmpi(subsmpl_size,'rand')
        sub_size = randi(data_dim_size,1); % random number
    else
        error('Wrong subsample size value!');
    end
else
    if subsmpl_size(1) < 1
        if length(subsmpl_size) == 1
            sub_size = round(subsmpl_size*data_dim_size);
        else
            a = subsmpl_size(1);
            b = subsmpl_size(2);
            randFromAB = a + (b-a).*rand();
            sub_size = round(randFromAB*data_dim_size);
        end
    else
        if length(subsmpl_size) == 1
            sub_size = subsmpl_size;
        else
            sub_size = randi([subsmpl_size(1),subsmpl_size(2)]);
        end
    end
end
end