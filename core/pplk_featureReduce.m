function [dataR,featInd] = pplk_featureReduce(data, method, k, varargin)
% [dataR,featInd] = pplk_featureReduce(data, method, k, varargin)
% Removes uninformative and redundant indices using feature reduction
% (unsupervised feature selection and extraction).
%
% INPUTS
%   data   
%       Matrix with unified validity indices.
%
%   method  
%       Method to use for feature reduction.
%
%           selection
%           - 'FSFS' 
%               Mitra et al., 2002.
%           - 'LS'   
%               Laplacian Score by He et al., 2005.
%           - 'SPEC' 
%               Spectral feature selection, Zhao & Liu 2007.
%           - 'FSKM'
%               Feature selection with k-medoids, Pepelka 2014.
%
%           extraction/transformation
%           - 'PCA'
%           - 'KPCA' 
%               Scholkopf et al., 1998.
%           - 'ICA'  
%               FastICA by Hyravinen et al., 2000.
%           - any of supported methods in the Dimensionality reduction
%             toolbox (by van der Maaten).
%             Unsupervised (prefered):
%               'ProbPCA', 'MDS', 'FactorAnalysis',  'Isomap', 'Laplacian',
%               'HessianLLE', 'LTSA','FastMVU', 'DiffusionMaps', 'SNE',
%               'SymSNE', 'tSNE', 'SPE', 'Autoencoder'
%             Unsupervised (occasional singularity problem on PRM data)
%               'GPLVM', 'Sammon', 'LandmarkIsomap', 'LLE', 'MVU', 'CCA',
%               'LandmarkMVU', 'LPP', 'NPE', 'LLTSA', 'LLC',
%               'ManifoldChart', 'CFA'
%             Supervised / labeled
%               'LDA', 'GDA', 'NCA', 'MCML', 'LMNN'
%           - 'FEKM' 
%               Feature extraction with k-means, Pepelka 2014.
%
%   k
%       - If k is a number, then it represents number of selected/extracted
%         features. Exception is the method 'FSFS', where k approximately
%         equals the number of selected features.
%
%       - If k is non-existent, empty or string, estimation of intrinsic
%         dimensionality is employed. If string, estimation is performed
%         using method with such a name.
%         Possible options are:
%           'CorrDim','NearNbDim','GMST','PackingNumbers','EigValue','MLE',
%           'MiND_ML','MiND_KL','DANCo','DANCoFit','kNN1','kNN2','kNN3',
%           'Hein','Takens'.
%               
%   varargin
%       Name-value list of optional parameters.
%           FSFS
%               'sim_method' 
%                   'cor' | 'lin' | 'mic'
%           LS
%               'neighbor_mode' 
%                   'KNN' | 'Supervised'
%               'weight_mode'   
%                   'Binary' | 'HeatKernel' | 'Cosine'
%               'kNN'           
%                   Number of k nearest neighbors.
%               'sigma'		
%                   If weight_mode is HeatKernel, sigma is used for
%                   gaussian kernel width.
%           KPCA
%               'kernel'
%                   'Gaussian' | 'Polynomial' | 'PolyPlus' | 'Linear'
%               'sigma'    
%                   Width of Gaussian kernel.
%               'd'    
%                   Param of Polynomial or PolyPlus kernel.
%           ICA
%               'kernel'	 
%                   'Gaussian' | 'pow3' | 'tanh' | 'skew'
%               'sigma'     
%                   Width of Gaussian kernel.
%               'd'      
%                   Param of tanh kernel.
%
%           OTHER:
%               varargin is directly passed to Dimensionality toolbox. See
%               implementation below for details.
%
%
% OUTPUTS
%   dataR
%       Reduced PRM matrix with the selected/extracted features.
%
%   featInd 
%       - Indices of the remaining features, if selection occures.
%       - Vector of NaN's (as many as number of remaining features) if
%         extraction occures.
%
%
% ACKNOWLEDGEMENTS AND REFERENCES
%
% Matlab Toolbox for Dimensionality Reduction
% (C)Laurens van der Maaten, Delft University of Technology
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you
% maintain the name of the original author.
%    L.J.P. van der Maaten, E.O. Postma, and H.J. van den Herik.
%      Dimensionality Reduction: A Comparative Review.
%      Tilburg University, Technical Report, TiCC-TR 2009-005, 2009.
%
% Feature Selection using Feature Similarity
%   P. Mitra, C. A. Murthy and S. K. Pal, Unsupervised Feature
%   Selection using Feature Similarity, IEEE Transactions on Pattern
%   Analysis and Machine Intelligence, Vol. 24, No. 4, pp 301-312,April 2002.
%   Code written by Pabitra Mitra.
%
% Laplacian Score
%   Xiaofei He, Deng Cai and Partha Niyogi, "Laplacian Score for
%   Feature Selection". Advances in Neural Information Processing
%   Systems 18 (NIPS 2005), Vancouver, Canada, 2005.
%   Code written by Deng Cai.
%
% Unsupervised Spectral feature selection
%   Zhao, Z., & Liu, H. (2007). Spectral feature selection for
%   supervised and unsupervised learning. Proceedings of the 24th
%   international conference on Machine learning (pp. 1151-1157)
%
% Kernel Principal Component Analysis
%   Bernhard Schölkopf, Alexander Smola, Klaus-Robert Müller. Nonlinear
%   Component Analysis as a Kernel Eigenvalue Problem", Neural Computation,
%   10:1299-1319, 1998.
%   Code written by Deng Cai (dengcai AT gmail.com).
%
% Independent Component Analysis
%   Hyvärinen, A., & Oja, E. (2000).
%   Independent component analysis: algorithms and applications.
%   Neural Networks, 13(4-5), 411-430.
%   FastICA code written by Hugo Gävert, Jarmo Hurri, Jaakko Särelä,
%   and Aapo Hyvärinen.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka
callDir=chdir(pplk_homeDir());

featInd = NaN;

dim = size(data,2);

% If not defined, estimate k with MLE
if ~exist('k','var') || isempty(k)
    k = pplk_dimEstimation(data,'MLE');
end
% If k is a string, use it as estimation method identifier
if ischar(k)
    k = pplk_dimEstimation(data,k);
end
if isnan(k) || isinf(k)
    warning('pplk:PRM','Unable to determine finite intrinsic dimensions. Setting k to size(data,2).');
    k=dim;
elseif k > dim
    warning('pplk:PRM','Estimated intrinsic dimensionality is larger than dimensionality of data. Setting k to size(data,2).');
    k=dim;
else
    k = round(k);
end
if k == 0
    warning('pplk:PRM','Estimated intrinsic dimensionality is ZERO. Setting k to size(data,2).');
    k = dim;
end

methodUp = upper(method);

if strcmpi(method,'NONE')
    dataR = data;
    featInd = 1:dim;
    chdir(callDir);
    return;
end

%fprintf(1,'Intrinsic dimensionality: %d\n',k);

switch methodUp
    
    % Feature Selection using Feature Similarity
    case 'FSFS'
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM']);
        
        if k == dim
            featInd = 1:k;
            dataR = data;
            return;
        end
        % parse varargin for parameters
        mask = find(strcmp('sim_method',varargin));
        if mask
            sim_method = varargin{mask+1};
            if ~any(ismember(sim_method,{'cor','lin','mic'}))
                warning('pplk_featureReduce:args','Wrong sim_method parameter! Using ''mic''.');
                sim_method = 'mic';
            end
        else
            sim_method = 'mic';
        end
        
        % s is scale parameter which decides the size of the reduced
        % feature set.
        % Approximately, s = #original features - k
        s = dim - k;
        [featInd, Rw] = FSFS(data,s,sim_method);
        % sort according to weights in Rw (descending)
        [~,indx]=sort(-Rw);
        featInd = featInd(indx);
        
%         if k ~= length(featInd)
%             fprintf(1,'\tFSFS: Output dimensions: %d\n',length(featInd));
%         end
        
        chdir(oldPath);
        dataR = data(:,featInd);
        
        % Laplacian Score - unsupervised feature selection
    case 'LS'
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'LS']);
        
        options = [];
        % parse varargin for parameters
        mask = find(strcmp('neighbor_mode',varargin));
        value = 'KNN';
        if mask
            value = varargin{mask+1};
            if ~any(ismember(value,{'KNN','Supervised'}))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        options.NeighborMode = value; % 'KNN' or 'Supervised'
        
        mask = find(strcmp('weight_mode',varargin));
        value = 'Cosine';
        if mask
            value = varargin{mask+1};
            if ~any(ismember(value,{'Binary','HeatKernel','Cosine'}))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        options.WeightMode = value; % 'Binary', 'HeatKernel' or 'Cosine'
        
        mask = find(strcmp('kNN',varargin));
        value = 5;
        if mask
            value = varargin{mask+1};
            if (value < 0) || (value > size(data,1))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        options.k = value;
        
        mask = find(strcmp('sigma',varargin));
        value = 1;
        if mask
            value = varargin{mask+1};
            if (value < 0) || (value > size(data,1))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        options.t = value;
        
        fea = NormalizeFea(data,2);
        W = constructW(fea,options);
        LS = LaplacianScore(data,W);
        
        [~, index] = sort(-LS); %larger the score, better the feature
        featInd = index(1:k);
        dataR = data(:,featInd);
        
        chdir(oldPath);
        
        % Unsupervised Spectral feature selection
    case 'SPEC'
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'SPEC']);
        
        fea = normData(data,2,0); % normalize instances, do not center data
        W = constructRBF(fea); % TODO: compare to constructW
        wFeat = fsSpectrum(W, data, 0);
        
        [~,index] = sort(wFeat); %smaller the score, better the feature
        featInd = index(1:k);
        dataR = data(:,featInd);
        chdir(oldPath);
        
        % Feature Selection using K-medoids (Nejc)
    case 'FSKM'
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'FKM']);        
        [dataR,featInd] = FSKM(data,k,50,0);        
        chdir(oldPath);
        
        % Feature Extraction using K-means (Nejc)
    case 'FEKM'
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'FKM']);        
        dataR = FEKM(data,k,50,0);
        featInd = nan(1,k);        
        chdir(oldPath);
        
        % Principal Component Analysis
    case 'PCA'
        % coeff stores eigenvectors
        % score is transformed data matrix
        % latent are eigenvalues
        [coeff,score,latent]=princomp(data);
        rankTolerance = 1e-7;
        nRedMax = sum (latent > rankTolerance);
        if nRedMax == 0,
            fprintf (['Eigenvalues of the covariance matrix are' ...
                ' all smaller than tolerance [ %g ].\n' ...
                'Please make sure that your data matrix contains' ...
                ' nonzero values.\nIf the values are very small,' ...
                ' try rescaling the data matrix.\n'], rankTolerance);
            error ('Unable to continue, aborting.');
        end
        if k>nRedMax
            warning('pplkPRM:PCA', ['Number of reduced features must be less than ',num2str(nRedMax),' (rank tolerance). Will return ', num2str(nRedMax), ' features.']);
            k = nRedMax;
        end
        dataR = score(:,1:k);
        featInd = nan(1,k);
        
        %         % Calculate PCA step by step
        %         % As done in FastICA, see also
        %         % http://matlabdatamining.blogspot.com/2010/02/principal-components-analysis.html
        %         % Center X by subtracting off column means
        %         PRM_c = bsxfun(@minus,PRM,mean(PRM,1));
        %         % Calculate the covariance matrix.
        %         PRM_cov = cov(PRM_c);
        %         % Calculate the eigenvalues and eigenvectors of covariance
        %         % matrix.
        %         [E, D] = eig(PRM_cov);
        %
        %         % The rank is determined from the eigenvalues - and not directly by
        %         % using the function rank - because function rank uses svd, which
        %         % in some cases gives a higher dimensionality than what can be used
        %         % with eig later on (eig then gives negative eigenvalues).
        %         rankTolerance = 1e-7;
        %         maxLastEig = sum (diag (D) > rankTolerance);
        %         if maxLastEig == 0,
        %             fprintf (['Eigenvalues of the covariance matrix are' ...
        %                 ' all smaller than tolerance [ %g ].\n' ...
        %                 'Please make sure that your data matrix contains' ...
        %                 ' nonzero values.\nIf the values are very small,' ...
        %                 ' try rescaling the data matrix.\n'], rankTolerance);
        %             error ('Unable to continue, aborting.');
        %         end
        %
        %         % Sort the eigenvalues - decending.
        %         eigenvalues = flipud(sort(diag(D)));
        %         data_pca = PRM_c * E(:,end:-1:1);
        
        
        % Kernel Principal Component Analysis
    case 'KPCA'
        nRedMax = dim;
        if k>nRedMax
            warning('pplkPRM:PCA', ['Number of reduced features must be less than ',num2str(nRedMax),'. Will return ', num2str(nRedMax), 'features.']);
            k = nRedMax;
        end
        
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM']);
        options = [];
        
        mask = find(strcmp('kernel',varargin));
        value = 'Gaussian';
        if mask
            value = varargin{mask+1};
            if ~any(ismember(value,{'Gaussian', 'Polynomial', 'PolyPlus', 'Linear'}))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        options.KernelType = value; % 'Gaussian' | 'Polynomial' | 'PolyPlus' | 'Linear'
        
        mask = find(strcmp('sigma',varargin));
        value = 1;
        if mask
            value = varargin{mask+1};
        end
        options.t = value;
        
        mask = find(strcmp('d',varargin));
        value = 2;
        if mask
            value = varargin{mask+1};
        end
        options.d = value;
        
        options.Kernel = 0;
        options.ReducedDim = k;
        
        [~,~,dataR] = kernelPCA(data,options);
        featInd = nan(1,k);
        
        chdir(oldPath);
        
        % Independent Components Analysis
    case 'ICA'
        % https://e-reports-ext.llnl.gov/pdf/240921.pdf , page 11:
        % "To find k < p independent components, one needs to first reduce the
        % dimension of the original data p to k, by a method such as PCA."
        oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'FastICA_25']);
        
        mask = find(strcmp('kernel',varargin));
        value = 'pow3';
        if mask
            value = varargin{mask+1};
            if ~any(ismember(value,{'Gaussian', 'pow3', 'tanh', 'skew'}))
                warning('pplk_featureReduce:args','Wrong parameter! Using defaults.');
            end
        end
        if strcmp(value,'Gaussian')
            value='gauss';
        end
        ICA_g = value; % 'gauss' | 'pow3' | 'tanh' | 'skew'
        
        mask = find(strcmp('sigma',varargin));
        value = 1;
        if mask
            value = varargin{mask+1};
        end
        ICA_a1 = 1/value;
        
        mask = find(strcmp('d',varargin));
        value = 1;
        if mask
            value = varargin{mask+1};
        end
        ICA_a2 = value;
        
        
        dataR = fastica(data', 'numOfIC', k, 'firstEig',1, 'lastEig',k, ...
            'approach','defl','g',ICA_g, 'a1',ICA_a1,'a2',ICA_a2, 'verbose', 'off');
        
        % if deflation approach did not converge, employ symmetric one
        if size(dataR,1) ~= k
            warning('pplk_featureReduce:ICA','ICA did not converge in deflation mode. Switching to symmetric.');
            dataR = fastica(data', 'numOfIC', k, 'firstEig',1, 'lastEig',k, ...
                'approach','symm','g',ICA_g, 'a1',ICA_a1,'a2',ICA_a2, 'verbose', 'off');
        end
        dataR = dataR';
        featInd = nan(1,k);
        chdir(oldPath);
        
    otherwise
        % include Dimension reduction toolbox (by van der Maaten)
        % Prefered methods are marked with *.
        % The parameters should be provided in this order, without names!
        %------------------------------------------------------------------
        %   METHOD          PARAMETERS (using default)
        %------------------------------------------------------------------
        %   *PCA:            - none
        %   *MDS:            - none
        %   *ProbPCA:        - <int> max_iterations -> default = 200
        %   *FactorAnalysis: - none
        %   GPLVM:          - <double> sigma -> default = 1.0
        %   Sammon:         - none
        %   *Isomap:         - <int> k -> default = 12
        %   LandmarkIsomap: - <int> k -> default = 12
        %                   - <double> percentage -> default = 0.2
        %   LLE:            - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   *Laplacian:      - <int> k -> default = 12
        %                   - <double> sigma -> default = 1.0
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   *HessianLLE:     - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   *LTSA:           - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   MVU:            - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   CCA:            - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   LandmarkMVU:    - <int> k -> default = 5
        %   *FastMVU:        - <int> k -> default = 5
        %                   - <logical> finetune -> default = true
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   *DiffusionMaps:  - <double> t -> default = 1.0
        %                   - <double> sigma -> default = 1.0
        %   *KernelPCA:      - <char[]> kernel -> {'linear', 'poly', ['gauss']}
        %                   - kernel parameters: type HELP GRAM for info
        %   *SNE:            - <double> perplexity -> default = 30
        %   *SymSNE:         - <double> perplexity -> default = 30
        %   *tSNE:           - <int> initial_dims -> default = 30
        %                   - <double> perplexity -> default = 30
        %   LPP:            - <int> k -> default = 12
        %                   - <double> sigma -> default = 1.0
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   NPE:            - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   LLTSA:          - <int> k -> default = 12
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   *SPE:            - <char[]> type -> {['Global'], 'Local'}
        %                   - if 'Local': <int> k -> default = 12
        %   *Autoencoder:    - <double> lambda -> default = 0
        %   LLC:            - <int> k -> default = 12
        %                   - <int> no_analyzers -> default = 20
        %                   - <int> max_iterations -> default = 200
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   ManifoldChart:  - <int> no_analyzers -> default = 40
        %                   - <int> max_iterations -> default = 200
        %                   - <char[]> eig_impl -> {['Matlab'], 'JDQR'}
        %   CFA:            - <int> no_analyzers -> default = 2
        %                   - <int> max_iterations -> default = 200
        
        %add Dimension reduction Toolbox to MATLAB path
        addpath(['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox'], ['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox',filesep,'techniques']);
        
        % use default parameters
        dataR = compute_mapping(data, method, k, varargin{:});
        featInd = nan(1,size(dataR,2));
        
        %remove Dimension reduction Toolbox from MATLAB path
        rmpath(['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox'], ['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox',filesep,'techniques']);
end

chdir(callDir);