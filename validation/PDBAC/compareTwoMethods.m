function [probGreater,diffMean,diffCI,C1,C2,diff] = compareTwoMethods(target,labels1,labels2,nsamples,alpha,res)
% Method with labels1 is better than method with labels2, if diffMean is
% positive (we can also consider confidence interval, so diffCI(1) has to be positive)
%
% probGreater - probability that method1 is better than method2

C1 = getcmClust(target,labels1);
C2 = getcmClust(target,labels2);


[probGreater,diff] = bac_dbp(C1,C2,nsamples,res);
diffMean = mean(diff);

diffCI = quantile(diff,[alpha/2 0.5 1-alpha/2]);
diffCI = [diffCI(1), diffCI(3)];
