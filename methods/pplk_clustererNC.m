function [labels, moreInfo]=pplk_clustererNC(data,K,params)
% VHODI:
%   data - data matrix; each row is an observation and each column a
%   feature
%   K - number of clusters
%   
%   Release notes:
%   2010, January 22: release of all c++ source mex files compatible with matlab R2009b
%   
%   Timothee Cour, Stella Yu, Jianbo Shi, 2004.

data=data';

% defaults
offset = 0.5;
scaleSigma = [];
if exist('params','var') && isstruct(params)
    if isfield(params,'NC_offset')
        offset = params.NC_offset;        
    end
    if isfield(params,'NC_scaleSigma')
        scaleSigma = params.NC_scaleSigma;        
    end
end

oldPath=chdir('NC');

%parameters
options = [];
options.scaleSigma = scaleSigma; % sigma for scaling the similarities between data
options.offset = offset; %offset in the diagonal of W, default 0.5
options.verbose = 0; %0 for verbose off mode, 1,2,3 for verbose on modes
options.maxiterations = 300; %max number of iterations in eigensolver
options.eigsErrorTolerance = 1e-8; %error tolerance in eigensolver
options.valeurMin=1e-6; %truncates any values in W less than valeurMin
% when computing similarity from distances, sigma is computed as: 
% scale_sig = 0.05*max(distances(:));

retVal=NC(data,K,options);

labels=retVal.target;
moreInfo.time=retVal.time;
moreInfo.options=retVal.options;

chdir(oldPath);
