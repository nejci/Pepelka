function [SD,SDbw]= indexSD_SDbw(data,labels,alpha)
% [SD,SDbw]= indexSD_SDbw(data,labels,alpha)
%--------------------------------------------------------------------------
% Cluster internal validity indices SD and S_Dbw (Halkidi et al.,2000&2001).
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   alpha       (scalar)    Required for SD.
%                           weighting factor (should equals Dis(labels(:,cmax))
%                           where cmax is the maximum number of input
%                           clusters). If empty set to number of clusters
%                           in labels.
%--------------------------------------------------------------------------
% OUTPUTS:
%   SD          (scalar)	value of the index SD
%   SDbw        (scalar)	value of the index SDbw
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Halkidi, M., Vazirgiannis, M., & Batistakis, Y. (2000). Quality scheme
% assessment in the clustering process. In Principles of Data Mining and
% Knowledge Discovery (pp. 265�276). London.
%
% Halkidi, M., & Vazirgiannis, M. (2001). Clustering validity assessment:
% finding the optimal partitioning of a data set. In Proceedings IEEE
% International Conference on Data Mining (pp. 187�194). Washington D.C.:
% IEEE Comput. Soc.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 28-May-2014 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[N,D] = size(data);
c = max(labels);

if ~exist('alpha','var') || isempty(alpha)
    alpha = c;
end


% compute the 2-norm of standard deviation of the data
stdev_S = sum(bsxfun(@minus,data,mean(data,1)).^2  ,1) ./ N;
stdev_S_norm = sqrt(stdev_S * stdev_S');

% compute centroids and standard deviation of clusters
v = zeros(c,D);
stdev_i_norm = zeros(c,1);
for i = 1:c
    C = data(labels==i,:);
    v(i,:) = mean(C,1);
    tmp = sum(bsxfun(@minus,C,v(i,:)).^2 ,1) ./ size(C,1);
    stdev_i_norm(i) = sqrt(tmp*tmp');
end

% compute average standard deviation for clusters
stdev_i_norm_sum = sum(stdev_i_norm);
stdev = sqrt(stdev_i_norm_sum) / c;

% compute average scattering for clusters - component Scat
Scat =  stdev_i_norm_sum / (c * stdev_S_norm);


% compute total separation between clusters Dis
dist_v = sqrt(sqdistance2(v,v));

Dmax = max(dist_v(:));
Dmin = min(dist_v(dist_v>0));
if isempty(Dmin)
    Dmin = 0;
end

Dis = 0;
for k=1:c    
    z=1:c;
    z(k)=[];
    Dis = Dis + 1 / sum(dist_v(k,z));
end
Dis = Dmax/Dmin * Dis;


% compute densities of custers
dens_v = zeros(c,1);
for i = 1:c    
    C = data(labels==i,:);
    dens_v(i) = density(C,v(i,:),stdev);
end

% compute inter-cluster densities - component Dens_bw
Dens = 0;
for i=1:c-1
    for j=i+1:c
        u_ij = (v(i,:)+v(j,:))/2;
        X = data(labels==i | labels == j,:);
        maxDens = max(dens_v(i),dens_v(j));
        if maxDens ~= 0
            Dens = Dens + density(X,u_ij,stdev) / maxDens;
        end        
    end
end
Dens_bw = 2*Dens/(c * (c-1));

SD = alpha * Scat + Dis;
SDbw = Scat + Dens_bw;


