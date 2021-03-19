% test, benchmark


N = 3000;
tol = 1e-10;

X = rand(N,1);
Y = rand(N,1);
data = [X,Y];

distMat = sqdistance2(data);


ticID = tic();
d1 = gabrielMatlab(distMat,tol);
t1 = toc(ticID);

ticID = tic();
d2 = gabriel(distMat,tol);
t2 = toc(ticID);

assert(isequal(d1,d2));
fprintf(1,'MATLAB: %f, MEX: %f\n',t1,t2);