% Higest Density Interval of the balanced accuracy.
%
% Usage:
%     [b_lower,b_upper] = bacc_hdi(C,alpha)
%
% Arguments:
% Arguments:
%     C - confusion matrix of classification outcomes
%     alpha - The HDI will cover 1-alpha of
%         probability mass such that (1-alpha)/2 remains on either end of
%         the distribution.
%
% Example:
%     C = [90, 10; 5 95];
%     alpha = 0.05;
%     [b_lower,b_upper] = bacc_hdi(C,alpha);
%
% Literature:
%
% [1] K. P. Murphy, Machine Learning: A Probabilistic Perspective (Adaptive
% Computation and Machine Learning series). Cambridge, MA:
% The MIT Press, 2012.
%
% Henry Carrillo, University of Zaragoza, Spain
% http://www.hcarrillo.com/
% -------------------------------------------------------------------------
function [PGreaterMC,diff] = bac_dbp(CA,CB,nsamples,res)

delta = .5;
proprnd = @(x) x + rand*2*delta - delta;

% Compute full distributions
AB = getParametersBac(CA);
pdf = @(x) betaavgpdf(x,AB,res);
x = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,'symmetric',1);

% Compute full distributions
AB1 = getParametersBac(CB);
pdf1 = @(x) betaavgpdf(x, AB1,res);
x1 = mhsample(1,nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);

diff = x-x1;
PGreaterMC = mean(x>x1);

end