function C = confMat3(target,predict)

n = length(target);
cat = spconvert([(1:n)' target ones(n,1)]);
cls = spconvert([(1:n)' predict ones(n,1)]);
cls = cls';
C = full(cls * cat);