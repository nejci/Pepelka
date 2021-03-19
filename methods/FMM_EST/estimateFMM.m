function [nComp,probs,means,covars] = estimateFMM(data,Kstart,method,params)
% numComp = onlineFMM(data,Kstart,alpha)
% Estimation of Finite Mixture Models from data
%--------------------------------------------------------------------------
% INPUTS		
%   data   (matrix)  input data, size of [N x D]
%   Kstart (scalar)  number of initial (max) finite mixture components
%                    default if emtpy: ceil(sqrt(N))
%   method (string)  'batch' - method of Figueiredo & Jain, 2002
%                    'online' - method of Zivkovic & Heijden, 2004
%   params (struct)  structure with optional fields:
%                      .replicates
%                        number of runs; return vector of estimations in
%                        nComp; probs, means and covars are cells
%
%                      .Kend (applies to batch)
%                        final number of components [default=1]
%
%                      .regularize (applies to batch)
%                        a regularization factor for covariance matrices; 
%                        in very small samples, it may be necessary to add 
%                        some small quantity to the diagonal of the covariances
%                        default if emtpy: 0
%
%                      .covarMode : (applies to batch)
%                        controls the covarince matrix:
%                          0: free covariances [default]
%                          1: diagonal covariances
%                          2: a common covariance for all components
%                          3: a common diagonal covarince for all components
%
%                      .alpha (applies to online)
%                             balances the influence of the data against
%                             the influence of the prior [default=1/N]
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
% Wrapper: Nejc Ilc, Pepelka toolbox
% References:
%   M. A. T. Figueiredo and A. K. Jain, "Unsupervised learning of
%   finite mixture models," IEEE Transactions on Pattern Analysis and
%   Machine Intelligence, vol. 24, no. 3, pp. 381–396, Mar. 2002.
%
%   Z. Zivkovic and F. van der Heijden, "Recursive unsupervised learning of
%   finite mixture models," IEEE transactions on pattern analysis and
%   machine intelligence, vol. 26, no. 5, pp. 651–6, May 2004.
%--------------------------------------------------------------------------

[N,d] = size(data);

% defaults
if ~exist('Kstart','var') || isempty(Kstart)
   Kstart = ceil(sqrt(N)); 
end

if ~exist('method','var') || isempty(method)
   method = 'batch'; 
end

replicates = 1;
Kend = 1;
regularize = 0;
covarMode = 0;
alpha = 1/N;

if exist('params','var') && isstruct(params)
    if isfield(params,'replicates') && ~isempty(params.replicates)
       replicates = params.replicates; 
    end
    if isfield(params,'Kend') && ~isempty(params.Kend)
       Kend = params.Kend; 
    end
    if isfield(params,'regularize') && ~isempty(params.regularize)
       regularize = params.regularize; 
    end
    if isfield(params,'covarMode') && ~isempty(params.covarMode)
       covarMode = params.covarMode; 
    end
    if isfield(params,'alpha') && ~isempty(params.alpha)
       alpha = params.alpha; 
    end
end



switch lower(method)
    case 'batch'
        parentDir = chdir('batch');        
        if replicates > 1
            nComp = zeros(1,replicates);
            probs = cell(1,replicates);
            means = cell(1,replicates);
            covars = cell(1,replicates);            
            for r = 1:replicates
                [nComp(r),probs{r},means{r},covars{r}] = batchFMM(data,Kstart,Kend,regularize,covarMode);
            end            
        else        
            [nComp,probs,means,covars] = batchFMM(data,Kstart,Kend,regularize,covarMode);
        end
        chdir(parentDir);
        
    case 'online'
        parentDir = chdir('online');
        if replicates > 1
            nComp = zeros(1,replicates);
            probs = cell(1,replicates);
            means = cell(1,replicates);
            covars = cell(1,replicates);            
            for r = 1:replicates
                [nComp(r),probs{r},means{r},covars{r}] = onlineFMM(data,Kstart,alpha);
            end            
        else        
            [nComp,probs,means,covars] = onlineFMM(data,Kstart,alpha);
        end
        chdir(parentDir);        
        
    otherwise
        error('Wrong estimation method!');        
end

