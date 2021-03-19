function [bmean,CI,chanceP,bmode,bmed,bnaive,bdist] = PDBAC(target,predict,alpha,show,res)
% Compute PDBAC (posterior distribution of the balanced accuracy) for
% confusion matrix C, given alpha.
%
% INPUTS:
%     target  - target class labels, a vector of values 1:g (g is number of
%               true classes); if predict is [], target is considered to be
%               the confusion matrix.
%     predict - predicted class labels by classification algorithm; a
%               vector of values 1:k (k is a number of predicted classes)
%     alpha   - the posterior probability interval will cover 1-alpha of
%               probability mass such that (1-alpha)/2 remains on either
%               end of the distribution.
%     show    - if 1 results are plotted.
%     res     - resolution/precision of sampling pdf
%
% OUTPUTS:
%     bmean   - Expected value of the balanced accuracy, i.e., the mean of
%               the average of two random variables that are independently
%               distributed according to Beta distributions.
%     CI      - two elements vector [low,high] of lower and upper bound of
%               confidence interval of the balanced accuracy.
%     chanceP - chance classifier threshold (1/n)
%     bmode   - Most likely balanced accuracy, i.e., the mode of the
%               average of 'n' random variables that are independently
%               distributed according to Beta distributions. Note that this
%               is NOT simply the mean of the modes of the individual
%               accuracies. This measure is bnaive.
%     bmed    - Median of the balanced accuracy, i.e., the median of the
%               average of 'n' random variables that are independently
%               distributed according to Beta distributions.
%     bnaive  - Naive balanced accuracy, i.e., simply the mean of the
%               individual accuracies - the mean of the accuracy modes.
%     bdist   - full distribution on the interval [0:res:1]
%
%--------------------------------------------------------------------------
% Writen by Nejc Ilc
% Based on the code from:
% Henry Carrillo (http://www.hcarrillo.com/) and
% Kay H. Brodersen (http://people.inf.ethz.ch/bkay/)
%--------------------------------------------------------------------------
% Reference:
% H. Carrillo, K. H. Brodersen, and J. A. Castellanos,
% "Probabilistic Performance Evaluation for Multiclass Classification Using
% the Posterior Balanced Accuracy," in in ROBOT2013: First Iberian Robotics
% Conference SE - 25, vol. 252, M. A. Armada, A. Sanfeliu, and M. Ferre,
% Eds. Springer International Publishing, 2014, pp. 347–361.

isConfusion = 0;
if ~exist('predict','var') || isempty(predict)
    isConfusion = 1;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('res','var') || isempty(res)
    res = 0.001;
end

% compute confusion matrix
if ~isConfusion
    C = getcm(target,predict);
else
    C = target;
end

n = size(C,1);
chanceP = 1/n;

AB = getParametersBac(C);
% Get means, modes, medians, PPIs, naive modes means of balanced accuracy (BAC)
% MEAN
%b = bacc_mean(C);
% Compute the sum using the characteristic function
% Set support
supportPdf = size(AB,1);
x = 0:res:supportPdf;
c = betaChfSum(res, AB);
bmean = sum(x.*c/supportPdf) * res;

if nargout > 1 || show
    % PPI - Posterior probability interval
    %[b_lower,b_upper] = bacc_ppi(C,alpha);
    b_lower = betaavginv(alpha/2,AB,res);
    b_upper = betaavginv(1-alpha/2,AB,res);
    CI = [b_lower,b_upper];
end

if nargout > 3 || show
    % MODE
    %bmode = bacc_mode(C);
    bmode = fminsearch(@(x) - betaavgpdf(x,AB,res), bmean);
end

if nargout > 4
    % MEDIAN
    %bmed = bacc_med(C);
    bmed = betaavginv(0.5, AB, res);
end

if nargout > 5
    % NAIVE - mean of accuracy modes
    %bnaive = bacc_naive(C);
    ABm1 = AB - 1;
    ABBac=0;
    % Set support
    supportPdf = size(ABm1,1);
    for i=1:1:supportPdf
        ABtemp = ABm1(i,1)/((ABm1(i,1)+ABm1(i,2)));
        if isnan(ABtemp)
            ABtemp = 0;
        end
        ABBac  = ABtemp+ABBac;
    end
    bnaive =(1/supportPdf)*ABBac;
end

if nargout > 6 || show
    % Compute full distributions
    x = 0:res:1;
    bdist = betaavgpdf(x, AB, res);
end

if show
    % plot it
    plotDistrBacIE2(AB,alpha,bmean,CI,chanceP,bmode,bdist,res,'northwest'); 
end