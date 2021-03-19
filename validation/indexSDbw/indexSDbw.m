function SDbw = indexSDbw(data,labels)
% SDbw = INDEXSDbw(data,labels)
%--------------------------------------------------------------------------
% Cluster internal validity index S_Dbw (Halkidi, 2001).
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%--------------------------------------------------------------------------
% OUTPUTS:
%   SDbw        (scalar)	value of the index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Halkidi, M., & Vazirgiannis, M. (2001). Clustering validity assessment:
% finding the optimal partitioning of a data set. In Proceedings IEEE
% International Conference on Data Mining (pp. 187�194). Washington D.C.:
% IEEE Comput. Soc. doi:10.1109/ICDM.2001.989517
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 15-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[N,D] = size(data);
[clusterIDs, ~, lbl] = unique(labels);
c = length(clusterIDs);


% compute the 2-norm of standard deviation of the data
stdev_S = sum(bsxfun(@minus,data,mean(data,1)).^2  ,1) ./ N;
stdev_S_norm = sqrt(stdev_S * stdev_S');

% compute centroids and standard deviation of clusters
v = zeros(c,D);
stdev_i_norm = zeros(c,1);
for i = 1:c
    C = data(lbl==i,:);
    v(i,:) = mean(C,1);
    tmp = sum(bsxfun(@minus,C,v(i,:)).^2 ,1) ./ size(C,1);
    stdev_i_norm(i) = sqrt(tmp*tmp');
end

% compute average standard deviation for clusters
stdev_i_norm_sum = sum(stdev_i_norm);
stdev = sqrt(stdev_i_norm_sum) / c;

% compute average scattering for clusters - component Scat
Scat =  stdev_i_norm_sum / (c * stdev_S_norm);

% compute densities of custers
dens_v = zeros(c,1);
for i = 1:c    
    C = data(lbl==i,:);
    dens_v(i) = density(C,v(i,:),stdev);
end

% compute inter-cluster densities - component Dens_bw
Dens = 0;
for i=1:c-1
    for j=i+1:c
        u_ij = (v(i,:)+v(j,:))/2;
        X = data(lbl==i | lbl == j,:);
        Dens = Dens + density(X,u_ij,stdev) / max(dens_v(i),dens_v(j) );
    end
end
Dens_bw = 2*Dens/(c * (c-1));

SDbw = Scat + Dens_bw;


