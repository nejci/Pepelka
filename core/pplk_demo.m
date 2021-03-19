function pplk_demo()
% Pepelka Package - Demoonstracija - Aleks

params = [];
[data, labelsT] = pplk_loadData(['artificial',filesep,'ring']);
[N,D] = size(data);
kT = length(unique(labelsT));

consMethod1 = 'STREHL-CSPA';
consMethod2 = 'LCE-CTS-SL';

options = [];
options.title = 'Ground truth of dataset Ring';
options.colorMode = 'color';
options.axisStyle = 'equal';
pplk_scatterPlot(data, labelsT, 2, options);

% {methodName, repetitions, k, k_mode}
k_sqrtN = floor(sqrt(N));
%methods = {'KM',2,[2, k_sqrtN],'rand'; ...
%		   'AL',2, [2, k_sqrtN], 'rand'};

clusterMethods = {
    'KM',8,[2, k_sqrtN],'rand';
    'SL', 4, [2, 5], 'rand'
    };

paramsGenEns = [];
paramsGenEns.subsampling = {'cols',[0.25,0.5]};
labelsEns = pplk_genEns(data, clusterMethods, paramsGenEns);

options.title = 'Ensemble members';
pplk_scatterPlot(data,labelsEns, [], options);

labelsCons1 = pplk_consEns(labelsEns, kT, consMethod1, params);
labelsCons2 = pplk_consEns(labelsEns, kT, consMethod2, params);

validInt1 = pplk_validInt(data, labelsCons1,{'SIL','DN','DNS'});
validExt1 = pplk_validExt(labelsT, labelsCons1, {'CA','AMI'});

validInt2 = pplk_validInt(data, labelsCons2,{'SIL','DN','DNS'});
validExt2 = pplk_validExt(labelsT, labelsCons2, {'CA','AMI'});

fprintf(1,'Cluster validity indices\n');
fprintf(1,'Index    \t%s\t%s\n', consMethod1, consMethod2);
fprintf(1,'Silhouete\t%f\t%f\n', validInt1.SIL, validInt2.SIL);
fprintf(1,'Dunn     \t%f\t%f\n', validInt1.DN, validInt2.DN);
fprintf(1,'Dunn mod \t%f\t%f\n', validInt1.DNS, validInt2.DNS);
fprintf(1,'-------------------------------------------------\n');
fprintf(1,'Clu. Acc.\t%f\t%f\n', validExt1.CA, validExt2.CA);
fprintf(1,'Adj. NMI \t%f\t%f\n', validExt1.AMI, validExt2.AMI);

options.title = ['Consensus func.: ', consMethod1, ', CA: ', num2str(validExt1.CA)];
pplk_scatterPlot(data, labelsCons1, [], options);

options.title = ['Consensus func.: ', consMethod2, ', CA: ', num2str(validExt2.CA)];
pplk_scatterPlot(data, labelsCons2, [], options);
