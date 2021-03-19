function [labels,moreInfo] = clusot(data,varargin)

% -clusot_surf-----------------------------------------------------------
% 'res' - (double) resolution of surface, default 0.05.
% 'type' - (string) computation method 'ellipse' or 'spline'. Default:
% 'ss_flag' - (bool) if set to true, special cases in spline
%             surface computation are handled by inserting mirrored
%             neuron positions. Default: true.
% 'flag' - (bool) if set to true, clusters are automerged if
%                 their cluster regions overlap in C. Default: false
% cluster_method - 'rec_flood' or 'clusot_cluster'
% Gradient based---------------------------------------------------------
% g - (double) threshhold parameter
% T - (double) resolution parameter for scanlines
% bflag - (bool) if set to true, clusters are automerged if
%                 their cluster regions overlap in C. Default: false.
% -rec_flood-------------------------------------------------------------
% 'theta' - (double) waterline threshold, applied w/r to current
%             local maximum, where 0 < theta <= 1. Default: 0.7.
% 'theta0' - (double) waterline threshold, which is used for
%             preflooding, where 0 <= theta0 < theta. Default: 0.3.
% 'tau' - (double) threshold for controlling the depth of 'valleys'
%             between clusters, where 0 <= tau < 1. Default: 0.

% Parameters defaults
% SOM
som_size = -1; % 'normal'
som_grid = 'hexa';
som_alg = 'batch';
som_trainlen = 'default';
% Clusot surface
type = 'spline';
res = 0.05;
ss_flag = true;
cluster_method = 'rec_flood';
flag = false;
% Gradient based
g = 0.5;
T = 256;
% Recursive flooding
theta0 = 0.3;
theta  = 0.7;
tau = 0;



% parse varargin options
i=1;
while i<=length(varargin),
    argok = 1;
    if ischar(varargin{i}),
        switch varargin{i},
            % argument IDs
            case 'som_size',  i=i+1; som_size = varargin{i};
            case 'som_grid',  i=i+1; som_grid = varargin{i};
            case 'som_alg',   i=i+1; som_alg = varargin{i};
            case 'som_trainlen',   i=i+1; som_trainlen = varargin{i};
            case 'res',       i=i+1; res = varargin{i};
            case 'type',      i=i+1; type = varargin{i};
            case 'ss_flag',   i=i+1; ss_flag = varargin{i};
            case 'cluster_method',   i=i+1; cluster_method = varargin{i};
            case 'g',         i=i+1; g = varargin{i};
            case 'T',         i=i+1; T = varargin{i};
            case 'flag',      i=i+1; flag = varargin{i};
            case 'theta',     i=i+1; theta = varargin{i};
            case 'theta0',    i=i+1; theta0 = varargin{i};
            case 'tau',       i=i+1; tau = varargin{i};
            otherwise
                argok=0;
        end
    else
        argok=0;
    end
    if ~argok,
        disp(['(clusot) Ignoring invalid argument #' num2str(i+1)]);
    end
    i = i+1;
end

% 1. level: SOM
path_somtoolbox = [pplk_homeDir('misc'),'somtoolbox_pplk'];
addpath(path_somtoolbox);

ticID = tic();
D = som_data_struct(data);
sM = som_make(D,...
    'msize',som_size,'lattice',som_grid,'algorithm',som_alg,'shape','sheet',...
    'neigh','gaussian','init','lininit','training',som_trainlen,'tracking',0);
tSOM = toc(ticID);

msize = sM.topol.msize;

% 2. level: Clusot
% If theta/theta0 or G are vectors, compute their values and validate with
% Davies-Bouldin index. Return result that maximize the index.

if strcmp(cluster_method,'rec_flood')
   [p,q] = meshgrid(theta0, theta);
   params = [p(:) q(:)];
   params = params(params(:,1) < params(:,2),:);
else
    params = g(:);
end

paramLen = size(params,1);
indexDB = nan(1,paramLen);
times = zeros(1,paramLen);
L = zeros(size(data,1),paramLen);
% loop over all parameters combinations
for p = 1:paramLen
    ticID = tic();
    if length(msize) == 2
        % run CLUSOT 2D
        if strcmp(cluster_method,'rec_flood')
%             clust = clusot2d(sM,D,'cluster_method',cluster_method,...
%                 'type',type,'res',res,'ss_flag',ss_flag,'flag',flag,...
%                 'theta0',params(p,1),'theta',params(p,2),'tau',tau);
            cmdStr = ['clust = clusotnd(sM,D,''type'',type,''res'',res,',...
                '''ss_flag'',ss_flag,''flag'',flag,''theta0'',params(p,1),',...
                '''theta'',params(p,2),''tau'',tau)'];
            evalc(cmdStr);
        else
            clust = clusot2d(sM,D,'cluster_method',cluster_method,...
            'type',type,'res',res,'ss_flag',ss_flag,'g',params(p),'T',T,'flag',flag);
            
        end
        
    else
        if ~strcmp(cluster_method,'rec_flood')
            error('Wrong cluster_method. Cannot apply to 3D SOM.');
        end
        % run CLUSOT 3D
        cmdStr = ['clust = clusotnd(sM,D,''type'',type,''res'',res,',...
            '''ss_flag'',ss_flag,''flag'',flag,''theta0'',params(p,1),',...
            '''theta'',params(p,2),''tau'',tau)'];
        evalc(cmdStr);
    end
    times(p) = toc(ticID);   
    labels = cluster_data(sM,D,clust);
    L(:,p) = labels;
    if paramLen > 1
        % compute Davies-Bouldin index
        if length(unique(labels)) > 1
            if min(labels) == 0
                labelsV = labels+1;
            else
                labelsV = labels;
            end
            [~,indexDB(p)] = pplk_validInt(D.data,labelsV,{'DB'});            
        end
    end
    
end
rmpath(path_somtoolbox);

[~,bestP] = min(indexDB);
labels = L(:,bestP);
tClusot = times(1);

moreInfo.numClust = max(labels);
moreInfo.time = tSOM+tClusot;
moreInfo.timeSOM = tSOM;
moreInfo.timeClusot = tClusot;
moreInfo.indexDB = indexDB;
moreInfo.paramsBest = params(bestP,:);
moreInfo.SOMstruct = sM;
moreInfo.params.som_size = som_size;
moreInfo.params.som_grid = som_grid;
moreInfo.params.som_alg = som_alg;
moreInfo.params.som_trainlen = som_trainlen;
moreInfo.params.type = type;
moreInfo.params.res = res;
moreInfo.params.ss_flag = ss_flag;
moreInfo.params.cluster_method = cluster_method;
moreInfo.params.flag = flag;
moreInfo.params.g = g;
moreInfo.params.T = T;
moreInfo.params.theta0 = theta0;
moreInfo.params.theta  = theta;
moreInfo.params.tau = tau;
