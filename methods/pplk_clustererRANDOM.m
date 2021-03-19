function [labels, moreInfo]=pplk_clustererRANDOM(data,K,params)

% [labels, moreInfo]=pplk_clustererRANDOM(data,K,params)
% returns random partition of the data into K clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modification: 7.3.2013 
% (C) Pepelka Package, Nejc Ilc (nejc.ilc@fri.uni-lj.si)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,D]=size(data);

tic 
labels=randi(K,N,1);
% ensure that there will be at least 1 sample with label from 1 to K
randInd = randsample(N,K);
labels(randInd) = 1:K;
time=toc;

moreInfo.time=time;

end