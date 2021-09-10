function dimEst = pplk_dimEstimation(data, method)
% dimEst = pplk_dimEstimation(data, method)
% Performs an estimation of the intrinsic dimensionality of data.
%
% INPUTS
%   data
%       A N-by-D matrix of data, where N is number of data samples and D is
%       number of dimensions.
%
%   method
%       Choose one of the following methods:
%
%       'CorrDim' 
%           Based on correlation dimension.
%       'NearNbDim' 
%           Based on nearest neighbor dimension.
%       'GMST' 
%           Based on the analysis of the geodesic minimum spanning tree.
%       'PackingNumbers' 
%           Based on the analysis of data packing number.
%       'EigValue' 
%           Based on analysis of PCA eigenvalues.
%       'MLE' 
%           Maximum likelihood estimator.
%       'MiND_ML' 
%          Minimum Neighbor Distance.
%       'MiND_KL' 
%           Minimum Neighbor Distance using Kullback-Leibler divergence.
%       'DANCo' 
%           Dimensionality from Angle and Norm Concentration.
%       'DANCoFit' 
%           Fast Dimensionality from Angle and Norm Concentration.
%       'kNN1', 'kNN2', 'kNN3' 
%           J. A. Costa and A. O Hero, 2003,2004.
%       'Hein' 
%          M. Hein & J-Y. Audibert, 2005.
%       'Takens'
%
%
% OUTPUTS
%   dimEst 
%       Estimated intrinsic dimensionallity of the data.
%
%
% ACKNOWLEDGEMENTS AND REFERENCES
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
% idEstimation toolbox
% URL: http://www.mathworks.com/matlabcentral/fileexchange/40112
% Copyright (c) 2013, Gabriele Lombardi. All rights reserved.
% Code covered by the BSD licence.
%   C.Ceruti, S.Bassis, A.Rozza, G.Lombardi, E.Casiraghi, P.Campadelli
%     DANCo: An intrinsic dimensionality estimator exploiting angle 
%     and norm concentration
%     Pattern Recognition, vol. 47, no. 8, pp. 2569-2581, 2014
%
% Matthias Hein and Jean-Yves Audibert
% URL: http://www.ml.uni-saarland.de/code/IntDim/IntDim.htm
%   M. Hein, J-Y. Audibert
%     Intrinsic dimensionality estimation of submanifolds in Euclidean
%     space Proceedings of the 22nd ICML, 289-296, Eds. L. de Raedt & S.
%     Wrobel,2005
%
% Alfred O. Hero
% URL: http://web.eecs.umich.edu/~hero/IntrinsicDim/
% Matab scripts for intrinsic dimension and entropy estimation using 
% k-nearest neighbor graphs. The details of the algorithms can be found in:
%    J. A. Costa and A. O Hero, "Entropic Graphs for Manifold Learning",
%      Proc. of IEEE Asilomar Conf. on Signals, Systems, and Computers,
%      Pacific Groove, CA, November, 2003.
%    J. A. Costa and A. O. Hero, "Geodesic Entropic Graphs for Dimension and
%      Entropy Estimation in Manifold Learning", to appear in IEEE Trans.
%      on Signal Processing, Aug., 2004.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka
callDir=chdir(pplk_homeDir());

if ~exist('method','var') || isempty(method)
    method = 'MLE';
end

% Dimensionality reduction toolbox
methodsDR = {'MLE','NearNbDim','CorrDim','PackingNumbers','GMST','EigValue'};
% idEstimation toolbox
methodsID = {'MiND_ML','MiND_KL','DANCo','DANCoFit'};
% Hein toolbox (needs compilation)
methodsHein = {'Hein', 'Takens'};
% Hero toolbox (needs compilation)
methodsHero = {'kNN1','kNN2','kNN3'};


% Preprocess data
% Remove duplicates from the dataset
% data = double(unique(data, 'rows'));
% 
% % Make sure data is zero mean, unit variance
% data = data - repmat(mean(data, 1), [size(data, 1) 1]);
% data = data ./ repmat(var(data, 1) + 1e-7, [size(data, 1) 1]);


if ismember(method, methodsDR)
    %add Dimension reduction Toolbox to MATLAB path
    addpath(['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox'], ['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox',filesep,'techniques']);
    
    dimEst = intrinsic_dim(data,method);
    
    %remove Dimension reduction Toolbox from MATLAB path
    rmpath(['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox'], ['..',filesep,'methods',filesep,'PRM',filesep,'drtoolbox',filesep,'techniques']);
    
elseif any(ismember(method,methodsID))
    oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'idEstimation']);
    
    funcName = str2func(method);
    dimEst = funcName(data');
    
    chdir(oldPath);
    
elseif any(ismember(method,methodsHein))
    oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'Hein']);
    
    out = GetDim(data');
    
    switch method
        case 'Hein'
            dimEst = out(1);
        case 'CorrDim'
            dimEst = out(2);
        case 'Takens'
            dimEst = out(3);
    end
    chdir(oldPath);
    
elseif any(ismember(method,methodsHero))
    oldPath = chdir(['..',filesep,'methods',filesep,'PRM',filesep,'kNN']);
    
    [n,d]=size(data); % number of sample points
    kneighbors=4; % number of neighbors in kNN
    % sample points to build bootstrap estimate of mean length function
    if n > (kneighbors+10)
        samp=(n-10):(n-1);
    else
        samp=(n-kneighbors+1):(n-1);
    end
    gamma=1;
    
    switch method
        case 'kNN1'
            % (M,N)=(1,10)
            dimEst = knn_graph_estim_1(data,kneighbors,gamma,1,10,samp);
        case 'kNN2'
            % (M,N)=(10,1)
            dimEst = knn_graph_estim_1(data,kneighbors,gamma,10,1,samp);
        case 'kNN3'
            % Dumb implementation of k-NN algorithm. Much slower in low
            % (extrinsic) dimensional spaces (like this example) but faster
            % for high (extrinsic) dimensional spaces (like images)
            
            % matrix of Euclidean distances between points
            Dist = L2_distance(data',data',1);
            % (M,N)=(1,10)
            dimEst = knn_graph_estim_2(Dist,kneighbors,gamma,1,10,samp);
            
    end
    chdir(oldPath);
    
else
    error('Wrong dimension estimation method!');
end

chdir(callDir);
