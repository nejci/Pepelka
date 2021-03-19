%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source code of the GP-MGLA method by Dong Huang.
% Version 1.0. May 14, 2014.
% Optimized for speed by Nejc Ilc, 10 June 2014
%
% If you use this code in your work, please cite the following paper:
%
% Dong Huang, Jian-Huang Lai, Chang-Dong Wang. 
% Combining Multiple Clusterings via Crowd Agreement Estimation and Multi-
% Granularity Link Analysis. Neurocomputing, in press, 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resultCls = gpmgla(baseCls, I_ncai, clsNums, alpha)

%% Build bipartite graph
%     disp('Build bipartite graph ... ');
[bcs, baseClsSegs] = getAllSegs(baseCls);
numClusters = length(baseClsSegs);
% disp(strcat('The total number of clusters in the ensemble:',num2str(numClusters)));
B = getBipartiteGraph(bcs, baseClsSegs, I_ncai, alpha); 

resultCls = [];
for i = 1:numel(clsNums); % clustering number. 
    if clsNums(i) > numClusters
        disp('Error in GP-MGLA:');
        disp('The cluster number in consensus clustering cannot exceeds the total number of clusters in the ensemble!!!');
        cls = ones(size(bcs,1),1); % set all labels to be one, indicating an error clustering.
    else
        cls = Tcut_for_clustering_ensemble(B,clsNums(i));
    end
    resultCls = [resultCls cls];
end 
clear labels B
