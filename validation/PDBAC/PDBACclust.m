function [PDBCA,CI,chanceP,bmode,bmed,bnaive,bdist,CA,C] = ...
    PDBACclust(target,predict,alpha,show,res)



if ~exist('predict','var') || isempty(predict)
    isConfusion = 1;
else
    isConfusion = 0;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('show','var') || isempty(show)
    show = 0;
end
if ~exist('res','var') || isempty(res)
    res = 0.001;
end

% compute confusion matrix
if ~isConfusion
    [C,~,cost] = getcmClust(target,predict);
else
    [C,~,cost] = getcmClust(target);
end

n = sum(C(:));
CA = -cost/n;

[PDBCA,CI,chanceP,bmode,bmed,bnaive,bdist] = PDBAC(C,[],alpha,show,res);
