function [labelsCons,numClust] = BCE(labelsEns,K,Ki,alphaInit,betaInit,lapParam)
% BCE - Bayesian Cluster Ensemble
% [labelsCons,numClust,weights] = WEAC(data,labelsEns,K,consFun,CVI,wMethod,CVIoptions)
% Ensemble clustering by Evidence accumulation clustering (Fred & Jain, 2005)
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns   (matrix)	clustering ensemble as [N X L] labels matrix;
%                           each column corresponds to one clustering. N is
%                           the number of data points and L the number of
%                           ensemble members.
%   K           (scalar)	number of clusters in output consensus
%	                        clustering; leave empty [] to automatically
%	                        determine this number by lifetime criterion.
%   Ki          (vector)    number of clusters of each ensemble member; can
%                           be left empty to compute it from labelsEns.
%   alphaInit   (vector)    vector [K X 1] of initial values for alpha. 
%                           If empty, they are randomly initialized.
%   betaInit    (cell)      cell [1 X L] of matrices with initial beta
%                           values. Each matrix is of size [K X Ki], where
%                           Ki is the number of clusters in i-th ensemble
%                           member. If empty, it is randomly initialized.
%   lapParam    (scalar)    parameter for Laplacian smoothing.
%                           
%
% OUTPUTS
%   labelsCons  (vector)    data labels - consensus clustering
%   numClust    (vector)    number of clusters in consensus clustering
%--------------------------------------------------------------------------
% EXAMPLE
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 3-October-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

[N,L] = size(labelsEns);

% number of clusters in ensemble members
if ~exist('Ki','var') || isempty(Ki)
    Ki = zeros(1,L);
    for iL = 1:L
        lbl = labelsEns(:,iL);
        Ki(iL) = length(unique(lbl));
    end
end

% parameter for laplace smoothing
if ~exist('lapParam','var') || isempty(lapParam)
    lapParam = 0.000001;
end

% Initial values for alpha
if ~exist('alphaInit','var') || isempty(alphaInit)
    alphaInit = rand(K,1);
end

% Initial values for beta
if ~exist('betaInit','var') || isempty(betaInit)
    betaInit = cell(1,L);
    for i=1:length(betaInit)
        temp = rand(K,Ki(i));
        temp = temp./(sum(temp,2)*ones(1,Ki(i)));
        betaInit{i} = temp;
    end
end


% learn BCE 
[phiAll,gammaAll,resultAlpha,resultBeta] = ...
    learnBCE(labelsEns,alphaInit,betaInit,lapParam,Ki);

% Obtain the cluster assignments from BCE
[~,labelsCons] = max(gammaAll',[],2);
numClust = max(labelsCons);



