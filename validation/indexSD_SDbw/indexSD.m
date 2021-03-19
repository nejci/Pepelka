function SD = indexSD(data,labels,alpha)
% SDbw = INDEXSD(data,labels,alpha)
%--------------------------------------------------------------------------
% Cluster internal validity index SD (Halkidi et al., 2000).
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   alpha       (scalar)    weighting factor (should equals to Dis(labels(:,cmax))
%                           where cmax is the maximum number of input
%                           clusters). If empty set to number of clusters
%                           in labels (not recommended).
%--------------------------------------------------------------------------
% OUTPUTS:
%   SD          (scalar)	value of the index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Halkidi, M., Vazirgiannis, M., & Batistakis, Y. (2000). Quality scheme
% assessment in the clustering process. In Principles of Data Mining and
% Knowledge Discovery (pp. 265�276). London.
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

% compute average scattering for clusters - component Scat
Scat =  sum(stdev_i_norm) / (c * stdev_S_norm);

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

SD = alpha * Scat + Dis;


