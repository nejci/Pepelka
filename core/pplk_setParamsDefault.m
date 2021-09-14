function params=pplk_setParamsDefault()
% params=pplk_setParamsDefault()
% Function returns default parameters' values for implemented clustering
% methods.
%
% INPUTS
%   none
%
%
% OUTPUTS		
%   params
%       Structure with default parameters' values.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka
%
% See also:		pplk_genEns

callDir=chdir(pplk_homeDir());

%% For all who care
params.distance = 'sqEuclidean'; % sqEuclidean, correlation
params.subsampling = {'none'};

%% AL, SL, WL
params.AL_distance='euclidean';
params.CL_distance='euclidean';
params.SL_distance='euclidean';
params.WL_distance='euclidean';

%% HCL
params.HCL_clustMethod = 'single';
params.HCL_distance = 'euclidean';

%% CLUSOT
params.CLUSOT_clustMethod = 'rec_flood';
params.CLUSOT_somSize = -1;
params.CLUSOT_somGrid = 'hexa';
params.CLUSOT_somAlg = 'batch';
params.CLUSOT_somTrainLen = 'default';
params.CLUSOT_res = 0.05;
params.CLUSOT_g = 0.5;
params.CLUSOT_theta0 = 0.3;
params.CLUSOT_theta  = 0.7;

%% CS
params.CS_sigma=[]; %Silverman's rule of thumb
params.CS_Kin=20;
params.CS_Nin=10;

%% EM GMM
clear eps;
params.EM_regularization = 1e-12; % A nonnegative regularization number added to the diagonal of covariance matrices to make them positive-definite.
params.EM_maxIter = 1000;

%% GSOM
params.GSOM_G = 0.001; %0.0008;
params.GSOM_dG = 0; %0.045;
params.GSOM_alfa = 0.01;
params.GSOM_maxIter = 0; % 0 = unlimited
params.GSOM_pOut = 0.1;
params.GSOM_msize = -2;%-1;
params.GSOM_shape = 'hexa';%'rect'; 
params.GSOM_showSOM = 0;
params.GSOM_showGrav = 0;
params.GSOM_distance='sqEuclidean';
params.GSOM_advanced={}; %additional settings for training parameters

%% SOMKM
params.SOMKM_nRuns = 5; 
params.SOMKM_msize = -2; % -2 means 2x the default size determined by somtoolbox
params.SOMKM_shape = 'rect';
params.SOMKM_run2toK = 0; % 0-run KM for K clusters, 1-run KM for range 2..K
                          % cluster and return solution which minimizes the
                          % Davies-Bouldin index
params.SOMKM_show = 0; % 1-visualize the results

%% SOMNC 
params.SOMNC_msize = -1; % default size determined by somtoolbox
params.SOMNC_shape = 'rect';

%% SOMSTAR
params.SOMSTAR_msize  = -1; % default size determined by somtoolbox
params.SOMSTAR_smooth = 2;  % amount of smoothing on U-matrix
params.SOMSTAR_show   = 0;  % show enhanced U-matrix and clustered data

%% KM
params.KM_maxIter=100; 
params.KM_nRuns=1;
params.KM_distance='sqEuclidean';

%% KVV NJW
params.KVV_sigma=0.2;
params.NJW_sigma=0.2;

%% NC
params.NC_scaleSigma = [];
params.NC_offset = 0.005;

%% SPECLS
params.SPECLS_Knn = 2; % in original paper it is 7
params.SPECLS_mode = 'LS';
params.SPECLS_KMruns = 20; % default 50
params.SPECLS_KMiters = 5; % default 10

%% Show results
params.showResult=0;

%% ENSEMBLE METHODS
% LCE
params.LCE_dc = 0.9; % decay factor; 0.9 is used in the original paper
params.LCE_R = 5; % number of iterations of SimRank algorithm
% COMUSA
% Relaxation of maximum edge weight constraint. 0.5 means 50% relaxation and
% thus bigger final clusters.
params.COMUSA_relaxation = 0.0; 
% HUANG
params.HUANG_betaOrINCAI = 2;
params.HUANG_alpha = 0.5;
% WEA
params.WEA_normalize = 1; % Should cluster properties be normalized?
% WEAC
params.WEAC_unifyMeth = 'range';
params.WEAC_reduceMeth = 'none';
params.WEAC_reduceDim = 'MLE';
params.WEAC_weightMeth = 'wMean2'; % aka. JWEAC
params.WEAC_weightMode = '';
% DICLENS-W
params.DICLENSW_unifyMeth = 'range';
params.DICLENSW_reduceMeth = 'none';
params.DICLENSW_reduceDim = 'MLE';
params.DICLENSW_weightMeth = 'wMean2'; % aka. JWEAC
params.DICLENSW_weightMode = '';

chdir(callDir);
