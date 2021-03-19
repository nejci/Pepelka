function [labels,moreInfo] = pplk_clustererSOMSTAR(data,K,params)
% Clustering of SOM with Starburst method (Hamel & Brown, 2011)
%--------------------------------------------------------------------------
% data - matrix with an observation in each row
% K - dummy; number of clusters is determined automatically
% params.SOMSTAR_
%    msize  - size of SOM grid; could be [], [xdim,ydim] or scaling factor -f
%    smooth - smoothing parameter (>0)
%    show   - if 1, plot the starburst map and clusters of data
%--------------------------------------------------------------------------
% Reference: 
% L. Hamel and C. Brown, "Improved Interpretability of the Unified Distance
% Matrix with Connected Components," in 7th International Conference on
% Data Mining, 2011, pp. 338–343.
%--------------------------------------------------------------------------
% Writen by Nejc Ilc
% Last modification 3-January-2014

% Defaults
somSize = -1;
smooth = [];
show = 0;
if exist('params','var') && isstruct(params)
    if isfield(params,'SOMSTAR_msize')
        somSize = params.SOMSTAR_msize;
    end
    if isfield(params,'SOMSTAR_smooth')
        smooth = params.SOMSTAR_smooth;
    end
    if isfield(params,'SOMSTAR_show')
        show = params.SOMSTAR_show;
    end
end

% run the SOMSTAR algorithm
algPath = pplk_homeDir('methods\SOMSTAR');
oldPath = chdir(algPath);

ticID = tic();
[labels,moreInfo] = SOMSTAR(data,somSize,smooth,show);
T = toc(ticID);

% remove path and cd back to parent dir
chdir(oldPath);

% return additional info
moreInfo.time = T;
