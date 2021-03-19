function [labels, moreInfo]=pplk_clustererFCM(data,K,params)
% Fuzzy c-means clustering 
% [center,U,obj_fcn] = fcm(data,cluster_n) 
% 
% [center, U, obj_fcn] = fcm(data, cluster_n) applies the fuzzy c-means
% clustering method to a given data set.
% 
% The input arguments of this function are:
%   data: data set to be clustered; each row is a sample data point
%   cluster_n: number of clusters (greater than one)
% 
% The output arguments of this function are
% 
%   center: matrix of final cluster centers where each row provides the
%           center coordinates
%   U: final fuzzy partition matrix (or membership function matrix)
%   obj_fcn: values of the objective function during iterations
% 
% fcm(data,cluster_n,options) uses an additional argument variable,
% options, to control clustering parameters, introduce a stopping criteria,
% set the iteration information display, or both.
% 
% options(1): exponent for the partition matrix U (default: 2.0)
% options(2): maximum number of iterations (default: 100)
% options(3): minimum amount of improvement (default: 1e-5)
% options(4): info display during iteration (default: 1)
% 
% If any entry of options is NaN, the default value for that option is used
% instead. The clustering process stops when the maximum number of
% iterations is reached or when the objective function improvement between
% two consecutive iterations is less than the minimum amount of improvement
% specified.

[N,D]=size(data);

tic

%default values - verbose off
options=[2.0, 100, 1e-5, 0];

[center,U,obj_fcn] = fcm(data,K,options);
maxU = max(U);

labels=zeros(N,1);
%index1 = find(U(1,:) == maxU);
%index2 = find(U(2, :) == maxU);
for ind=1:K
    labels(U(ind,:) == maxU)=ind;
end
time=toc;

moreInfo.time=time;