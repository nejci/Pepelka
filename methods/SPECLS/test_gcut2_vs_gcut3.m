% benchmark gcut2 vs. gcut3

KMruns = 20;
KMiters = 5;

load precomputedV_knn2.mat
data_num = length(datasets);

T_true = zeros(data_num,2);

S_true = zeros(data_num,2);

for data_i = 1:data_num
    data_el = datasets{data_i};
    fprintf(1,'%s\n',data_el);
    % load precomputed SVD
    Vd = V{data_i};
    % load data and ground-truth labels
    [data,target] = pplk_loadData(data_el);
    KTRUE = max(target);
    
    ticID = tic();
    labels = gcut2(Vd,KTRUE,KMruns,KMiters);
    T_true(data_i,1) = toc(ticID);
    [~,S_true(data_i,1)] = pplk_validExt(target,labels,{'CA'});
    fprintf(1,'\t%f s (%f)\n',T_true(data_i,1),S_true(data_i,1));
    
    ticID = tic();
    labels = gcut3(Vd,KTRUE,KMruns,KMiters);
    T_true(data_i,2) = toc(ticID);    
    [~,S_true(data_i,2)] = pplk_validExt(target,labels,{'CA'});
    fprintf(1,'\t%f s (%f)\n',T_true(data_i,2),S_true(data_i,2));
        
end

RES = [datasets, num2cell(T_true), num2cell(S_true), num2cell(S_true(:,1) - S_true(:,2))]
[sum(T_true),mean(S_true)]