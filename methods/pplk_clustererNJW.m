function [labels, moreInfo]=pplk_clustererNJW(data,K,params)
%nisem zadovoljen z rezultati, sigma ni ok
% VHODI:
%   data - matrika vhodnih podatkov (velikosti Nxd, kjer je N število vzorcev in d njihova razsežnost)
%   K - število gruè, v katere razvršèamo vzorce 
%

% defaults
sigma = 0.2;
if exist('params','var') && isstruct(params)
    if isfield(params,'NJW_sigma')
        sigma = params.NJW_sigma;        
    end
end

oldPath=chdir(['..',filesep,'libs',filesep,'SpectraLib']);

%init toolbox (path to SpectraLib is added (SPECTRAL_HOME) and some global
%variables are created)
startup


% S_ij= exp(-distance_ij^2/(2*sigma*sigma));
S=S_from_points(data',sigma,0,0);

tic
labels=njw_kmeans(S,K);
time=toc;

moreInfo.time=time;
moreInfo.options.sigma=sigma;

%clear path and global variables created by startup script
rmpath(ADDED_PATH);
clearGlobals;
chdir(oldPath);