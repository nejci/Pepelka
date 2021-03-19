function [labels, moreInfo]=pplk_clustererKVV(data,K,params)
% VHODI:
%   data - matrika vhodnih podatkov (velikosti Nxd, kjer je N �tevilo vzorcev in d njihova razse�nost)
%   K - �tevilo gru�, v katere razvr��amo vzorce 
%

% defaults
sigma = 0.2;

if exist('params','var') && isstruct(params)
    if isfield(params,'KVV_sigma')
        sigma = params.KVV_sigma;        
    end
end

oldPath=chdir('..\misc\SpectraLib');

%init toolbox (path to SpectraLib is added (SPECTRAL_HOME) and some global variables are created)
startup;

% S_ij= exp(-distance_ij^2/(2*sigma*sigma));
S = S_from_points(data',sigma,0,0);

tic
labels=kvv1_ncut(S,K);
time=toc;

moreInfo.time=time;
moreInfo.sigma=sigma;

%clear path and global variables created by startup script
rmpath(ADDED_PATH);
clearGlobals;
chdir(oldPath);
