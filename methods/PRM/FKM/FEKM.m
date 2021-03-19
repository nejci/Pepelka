function [dataR,labels,energy] = FEKM(data,k,iters,normalize)
% Feature Extraction using K-Means

if ~exist('iters','var') || isempty(iters)
    iters = 50;
end
if ~exist('normalize','var') || isempty(normalize)
    normalize = 0;
elseif normalize == 1
    normalize = 'range';    
end

feats = data';
if normalize
    feats = pplk_normalize(feats,normalize);
end

[labels,C,sumd,D] = kmeans(feats,k,...
    'distance','sqEuclidean',...
    'emptyaction','drop',...
    'replicates',iters,...
    'start','sample');
energy = sum(sumd);
%fprintf(1,'k=%d, ssumd: %f\n',k,energy);
dataR = C';

% figure();
% plot(feats','k-.');
% hold on
% plot(dataR,'r');
% 
% figure();
% pplk_scatterPlot(feats,labels);