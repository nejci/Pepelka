function [dataR,featInd,label,energy] = FSKM(data,k,iters,normalize)
% Feature Selection using K-Medoids

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

energy = Inf;
label = [];
featInd = [];

% repeat iters-times and pick the solution which minimizes energy
for r = 1:iters
    [label_i, energy_i, featInd_i] = kmedoids(feats',k);
    if energy_i < energy
        label = label_i;
        energy = energy_i;
        featInd = featInd_i;
    end
end

%fprintf(1,'k=%d, energy: %f\n',k,energy);
dataR = feats(featInd,:)';

% figure();
% plot(feats','k-.');
% hold on
% plot(dataR,'r');
% 
% figure();
% pplk_scatterPlot(feats,label);