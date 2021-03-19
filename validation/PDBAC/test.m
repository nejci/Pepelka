
%% build test data
K = 2:20;
N = [20,30,50,100,200,500,1000];
nK = length(K);
nN = length(N);
TARGET = cell(nK,nN);
PREDICT0 = cell(nK,nN);
PREDICT1 = cell(nK,nN);
PREDICT2 = cell(nK,nN);

ki = 1;
for k = K
    ni = 1;
    for n = N
        target = ones(n,1);
        bin = floor(n/k);
        for c = 2:k
            target(bin*(c-1)+1:bin*c,1) = c;
        end
        
        %predict = randi(k,[n,1]);
        predictT = target;
        predict = predictT;
        perm = randperm(k);
        for c = 1:k
            predictT(predict==c) = perm(c);
        end
        
        TARGET{ki,ni} = target;
        PREDICT0{ki,ni} = ones(n,1);
        PREDICT1{ki,ni} = predictT;
        PREDICT2{ki,ni} = randi(k,[n,1]);
        
        ni = ni+1;
    end
    ki = ki+1; 
end

%% test
[nK,nN] = size(TARGET);
alpha = 0.05;
for ki = 1:nK
    for ni = 1:nN
        fprintf(1,'ki=%d, ni=%d\n',ki,ni);
        tic();
        [BCA1_0, CI,chanceP,bmode,bmed,bnaive,bdist,CA,C0] = PDBACclust(TARGET{ki,ni},PREDICT0{ki,ni},alpha,0,0.001);
        [BCA1_1, CI,chanceP,bmode,bmed,bnaive,bdist,CA,C1] = PDBACclust(TARGET{ki,ni},PREDICT1{ki,ni},alpha,0,0.001);
        [BCA1_2, CI,chanceP,bmode,bmed,bnaive,bdist,CA,C2] = PDBACclust(TARGET{ki,ni},PREDICT2{ki,ni},alpha,0,0.001);
        t1 = toc();
        fprintf(1,'    1. OK! time: %f\n',t1);
        
        oldpath = chdir('..\core\bac');
        tic();
        [BCA2_0,CI,chanceP,bmode,bmed,bnaive,y_bacc] = PDBAorig(C0, alpha);
        [BCA2_1,CI,chanceP,bmode,bmed,bnaive,y_bacc] = PDBAorig(C1, alpha);
        [BCA2_2,CI,chanceP,bmode,bmed,bnaive,y_bacc] = PDBAorig(C2, alpha);
        t2 = toc();
        chdir(oldpath);
        fprintf(1,'    time2: %f\n',t2);
        
        assert(BCA1_0 == BCA2_0, ['0: ki=',num2str(ki),', ni: ', num2str(ni)]);
        assert(BCA1_1 == BCA2_1, ['1: ki=',num2str(ki),', ni: ', num2str(ni)]);
        assert(BCA1_2 == BCA2_2, ['2: ki=',num2str(ki),', ni: ', num2str(ni)]);
        fprintf(1,'    test passed!\n');
        
        
    end
end

disp('Test successfully passed!');