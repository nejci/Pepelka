% test EM

datasets = {'iris'};

for d = 1:length(datasets)
    dataset_el = datasets{d};
    
    [data, target] = pplk_loadData(dataset_el);
    data = pplk_normalize(data,'zscore');
    
    K = max(target);
    
    regularization = 1e-12;
    maxIter = 1000;
    
    % MATLAB Statistics toolbox
    ticID = tic();
    labels1 = EM_statToolbox(data,K,regularization,maxIter);
    t1 = toc(ticID);
    %pplk_scatterPlot(data,labels1);
    
    % Mo Chen implementation
    ticID = tic();
    labels2 = emgm(data',K,regularization,maxIter);
    t2 = toc(ticID);
    %pplk_scatterPlot(data,labels2);
    
    fprintf(1,'%s: %f | %f\n', dataset_el, t1, t2);
end