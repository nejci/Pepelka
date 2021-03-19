% convert model file
% function handles have also data stored in them!!!

A = load('DANCo_fits.mat');

f1 = functions(A.fitDhat);
f2 = functions(A.fitMu);
f3 = functions(A.fitTau);

% extract data
fitDhatFun = f1.workspace{1}.fitDhatFun;
fitMuFun = f2.workspace{1}.fitMuFun;
fitTauFun = f3.workspace{1}.fitTauFun;
% create new anonymous functions
fitDhat = @(D,N)(fnval2(fitDhatFun,{D,N})');
fitMu = @(D,N)(fnval2(fitMuFun,{D,N})');
fitTau = @(D,N)(fnval2(fitTauFun,{D,N})');
k = A.k;

% Save to new file
modelfile = 'DANCo_fits2.mat';
save(modelfile,'fitMu','fitTau','fitDhat','k');