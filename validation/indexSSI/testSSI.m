% test of SSI with kmeans on Iris

% [dataOrig, target] = pplk_loadData('real\iris_orig');
% data = dataOrig;
% 
% % cluster data
% res = [];
% resW = [];
% Kmax = 12;
% for K = 2:Kmax
%     labels = pplk_runClusterer('KM',data,K,1);
% 
%     [SSI, SSIW] = indexSSI(data,labels,0);
%     res = [res, SSI];
%     resW = [resW, SSIW];
%     
% end
% [res;resW]

% first run runSSI_cclust.R from
% D:\Nejc\LASPP\mojeDelo\Pepelka\validation\test_compare_with_R\
% It will create two files with centers and labels. Use them to validate
% index.
centers = load('D:\Nejc\LASPP\mojeDelo\Pepelka\validation\test_compare_with_R\SSI_centers.txt');
labels = load('D:\Nejc\LASPP\mojeDelo\Pepelka\validation\test_compare_with_R\SSI_labels.txt');
[SSI, SSIW] = indexSSI(centers,labels,1)