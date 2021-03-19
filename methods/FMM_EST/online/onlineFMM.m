function [nComp,probs,means,covars] = onlineFMM(data,Kstart,alpha)
% [numComp,probs,means,covars] = onlineFMM(data,Kstart,alpha)
% Recursive Unsupervised Learning of Finite Mixture Models (Zivkovic &
% Heijden, 2004)
%--------------------------------------------------------------------------
% INPUTS		
%   data   (matrix)  input data, size of [N x D]
%   Kstart (scalar)  number of initial finite mixture components 
%                    default if emtpy: ceil(sqrt(N))
%   alpha  (scalar)  parameter that balances the influence of the data
%                    against the influence of the prior;
%                    default if emtpy: 1/N
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
% Author:   Z.Zivkovic (10-2-2003), http://www.zoranz.net/
% Wrapper: Nejc Ilc, Pepelka toolbox
% Reference:
%   Z. Zivkovic and F. van der Heijden, "Recursive unsupervised learning of
%   finite mixture models," IEEE transactions on pattern analysis and
%   machine intelligence, vol. 26, no. 5, pp. 651–6, May 2004.
%--------------------------------------------------------------------------
data = data';
N = size(data,2);

if ~exist('Kstart','var') || isempty(Kstart)
   Kstart = ceil(sqrt(N)); 
end

if ~exist('alpha','var') || isempty(alpha)
   alpha = 1/N; 
end


[probs,means,covars] = EMRandomInit(data,Kstart);

for iData=1:N
    %algorithm - update step
    [probs,means,covars] = EMStep(probs,means,covars,data(:,iData),alpha);
    
end

nComp = sum(probs>0);
