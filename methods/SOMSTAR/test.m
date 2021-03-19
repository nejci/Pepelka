% Test of SOMSTAR

% load data
[data,target] = pplk_loadData('iris');
data = pplk_normalize(data,'zscore');

% run SOMSTAR directly
somSize = [15,10];
smooth = 0;
show = 1;
[labels,moreInfo] = SOMSTAR(data,somSize,smooth,show);

% run wrapper
params = pplk_setParamsDefault();
params.SOMSTAR_msize = somSize;
params.SOMSTAR_smooth = smooth;
params.SOMSTAR_show = 1;
[labels2,moreInfo2] = pplk_runClusterer('SOMSTAR',data,[],1,params);