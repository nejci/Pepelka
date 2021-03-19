



N = 5;
K = 2;
M = 3;
E = randi(K,N,M);
%E = [1 1;2 2;1 2];
dc = 0.8;


[Enew, no_allcl, C] = relabelCl(E); % re-labelling clusters in the ensemble E

tic;
[wcl0 pc0]= weightCl0(Enew);
toc

mex -largeArrayDims OPTIMFLAGS="/openmp $OPTIMFLAGS" weightCl.c
tic
[wcl_mex pc_mex]= weightCl_mex(Enew,no_allcl);
toc
assert(isequal(wcl0,wcl_mex),'wcl differ!');
assert(isequal(pc0,pc_mex),'pc differ!');

%% -------------------------------------------------------------------------

N = 1000;
K = 20;
M = 20;
E = randi(K,N,M);
dc = 0.8;

%%
tic
S0 = cts0(E, dc);
toc

%% 
mex -largeArrayDims cts_S_mex.c OPTIMFLAGS="/openmp $OPTIMFLAGS"

iter = 10;
t_mex = zeros(1,iter);
for i = 1:iter    
    tic
    S_mex = cts(E,dc);
    t_mex(i) = toc();
    
    diff = sum(sum(abs(S0-S_mex)));
    if  diff > 1e-9
        error('S differ!');
    end 
end
mean(t_mex)