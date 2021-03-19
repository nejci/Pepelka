%% -------------------------------------------------------------------------

N = 1000;
K = 20;
M = 20;
E = randi(K,N,M);
%E = [1 1 1;1 2 2;2 2 3;2 3 4];
% Toy example from paper
%E = [1 1 2 2 3 3; 3 3 2 2 1 1; 1 1 1 2 2 2; 1 1 nan nan nan 2]';

%%

ticID = tic();
[labelsCons2, Kcons2] = diclens(E);
t2 = toc(ticID);

oldDir = chdir('..\fast_old');
ticID = tic();
[labelsCons1, Kcons1] = diclens(E);
t1 = toc(ticID);
chdir(oldDir);



fprintf(1,'Fast:  %f\nFast2: %f\n---------\n',t1,t2);
assert(isequal(labelsCons1,labelsCons2),'!!!');
assert(isequal(Kcons1,Kcons2),'!!!');
assert(isequaln(avgICS1,avgICS2),'!!!');
assert(isequaln(avgECS1,avgECS2),'!!!');
assert(isequal(finalClusters_hist1,finalClusters_hist2),'!!!');