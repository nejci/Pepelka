% test of Huang functions

%% load data
[data target] = pplk_loadData('iris');
[N D] = size(data);
K = max(target);

%% Generate ensemble
methods = {'KM', 10, [2 2*floor(sqrt(N))], 'rand'};
params = pplk_setParamsDefault();
labelsEns = pplk_genEns(data,methods,params);

%% Consensus
betaVal = 2;
alphaVal = [];

ticID = tic();
ncai = getNCAI(labelsEns);
% The influence of the NCAI
I_ncai = ncai.^betaVal;

params.HUANG_betaOrINCAI = I_ncai;
params.HUANG_alpha = alphaVal;

[labelsCons_WEAC_SL, Kcons_WEAC_SL] = pplk_consEns(labelsEns, K, 'HUANG-SL', params);
[labelsCons_WEAC_CL, Kcons_WEAC_CL] = pplk_consEns(labelsEns, K, 'HUANG-CL', params);
[labelsCons_WEAC_AL, Kcons_WEAC_AL] = pplk_consEns(labelsEns, K, 'HUANG-AL', params);
[labelsCons_MGLA_GP, Kcons_MGLA_GP] = pplk_consEns(labelsEns, K, 'HUANG-GPMGLA', params);

time = toc(ticID);

%% Validate
validExt_WEAC_SL = pplk_validExt(target, labelsCons_WEAC_SL, {'NMI','CA','BCA'});
validExt_WEAC_CL = pplk_validExt(target, labelsCons_WEAC_CL, {'NMI','CA','BCA'});
validExt_WEAC_AL = pplk_validExt(target, labelsCons_WEAC_AL, {'NMI','CA','BCA'});
validExt_MGLA_GP= pplk_validExt(target, labelsCons_MGLA_GP, {'NMI','CA','BCA'});

fprintf(1,'Done! Time elapsed: %f\n',time);
fprintf(1,'NMI scores\n');
fprintf(1,'WEAC-SL: %f\n',validExt_WEAC_SL.NMI);
fprintf(1,'WEAC-CL: %f\n',validExt_WEAC_CL.NMI);
fprintf(1,'WEAC-AL: %f\n',validExt_WEAC_AL.NMI);
fprintf(1,'GP-MGLA: %f\n',validExt_MGLA_GP.NMI);
