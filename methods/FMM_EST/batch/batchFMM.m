function [nComp,probs,means,covars] = batchFMM(data,Kstart,Kend,regFact,covarMode)
% [numComp,probs,means,covars] = batchFMM(data,Kstart,Kend,regularize,covarMode)
% Recursive Unsupervised Learning of Finite Mixture Models (Zivkovic &
% Heijden, 2004)
%--------------------------------------------------------------------------
% INPUTS		
%   data   (matrix)  input data, size of [N x d]
%   Kstart (scalar)  initial (max) number of components to be tested
%                    default if emtpy: ceil(sqrt(N))
%   Kend   (scalar)  final (min) number of components to be tested
%                    default if emtpy: 1
%   regFact(scalar)  a regularization factor for covariance matrices; in
%                    very small samples, it may be necessary to add some
%                    small quantity to the diagonal of the covariances
%                    default if emtpy: 0
%
%   covarMode        controls the covarince matrix
%                      0: free covariances
%                      1: diagonal covariances
%                      2: a common covariance for all components
%                      3: a common diagonal covarince for all components
%
%
% OUTPUTS
%   nComp  (scalar)  number of mixture components
%   probs  (vector)  vector of mixture probabilities
%   means  (matrix)  the estimates of the means of the components [d X nComp]
%   covars (array)   contains the estimates of the covariances of the components
%                    It is 3D array, such that covars(:,:,1) is the d X d
%                    covariance of the first component and covars(:,:,nComp) 
%                    is the covariance of the last component
%--------------------------------------------------------------------------
% Author:  Mario Figueiredo (mtf@lx.it.pt), http://www.lx.it.pt/~mtf/
% Wrapper: Nejc Ilc, Pepelka toolbox
% Reference:
%   M. A. T. Figueiredo and A. K. Jain, "Unsupervised learning of
%   finite mixture models," IEEE Transactions on Pattern Analysis and
%   Machine Intelligence, vol. 24, no. 3, pp. 381–396, Mar. 2002.
%--------------------------------------------------------------------------
data = data';
N = size(data,2);

if ~exist('covarMode','var') || isempty(covarMode)
   covarMode = 0; 
end
if ~exist('regularize','var') || isempty(regFact)
   regFact = 0; 
end
if ~exist('Kend','var') || isempty(Kend)
   Kend = 1; 
end
if ~exist('Kstart','var') || isempty(Kstart)
   Kstart = ceil(sqrt(N)); 
end

% stopping threshold
threshold = 1e-4;


[nComp,probs,means,covars,dl,countf] = ...
    mixtures4(data,Kend,Kstart,regFact,threshold,covarMode);
