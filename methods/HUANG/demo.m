%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The demo of the WEAC and GP-MGLA methods by Dong Huang.
% Version 1.0. May 14, 2014.
%
% If you use this code in your work, please cite the following paper:
%
% Dong Huang, Jian-Huang Lai, Chang-Dong Wang. 
% Combining Multiple Clusterings via Crowd Agreement Estimation and Multi-
% Granularity Link Analysis. Neurocomputing, in press, 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters  
para.alpha = 0.5;
para.beta = 2;

%% Load the base clusterings and the groundtruth
% baseCls is a n x m matrix. Each column in baseCls represents a base
% clustering. There are totally m base clusterings.
load('baseCls_and_gt.mat','baseCls','gt');

%% Get the NCAI
ncai = getNCAI(baseCls);
% The influence of the NCAI
I_ncai = ncai.^para.beta;

%% Seting number of clusters
% clsNums is a 1xd vector. E.g., with clsNums=[2,3,4,5,9,20], the weac method 
% will output 6 consensus clusterings that contains 2, 3, 4, 5, 9, and 20 
% clusters respectively.
clsNums = [2, 3, 4, 5, 9, 20];

%% Perform WEAC and GP-MGLA
% The resultCls_sl is the consensus result of WEAC-SL, while resultCls_cl
% WEAC-CL and resultCls_al WEAC-AL.
[resultCls_sl, resultCls_cl, resultCls_al] = weac(baseCls, I_ncai, clsNums); 

% Compute the consensus clustering by GP-MGLA.
resultCls_gp = gpmgla(baseCls, I_ncai, clsNums, para.alpha);

% Pepelka interface methods
% Kind = 6;
% [labelsCons1 numClust1] = HUANG_WEAC(baseCls, clsNums(Kind), 'single', I_ncai);
% [labelsCons2 numClust2] = HUANG_GPMGLA(baseCls, clsNums(Kind), para.beta, para.alpha);
% 
% isequal(resultCls_sl(:,Kind), labelsCons1)
% isequal(resultCls_gp(:,Kind), labelsCons2)

%% Display the NMI scores of consensus results.
disp('The NMI of the consensus results of WEAC-AL:');
for i = 1:size(resultCls_al,2)
    nmi = NMImax(resultCls_al(:,i),gt);
    disp(['k=',num2str(clsNums(i)),' : NMI=', num2str(nmi)]);
end

disp('The NMI of the consensus results of GP-MGLA:');
for i = 1:size(resultCls_al,2)
    nmi = NMImax(resultCls_gp(:,i),gt);
    disp(['k=',num2str(clsNums(i)),' : NMI=', num2str(nmi)]);
end