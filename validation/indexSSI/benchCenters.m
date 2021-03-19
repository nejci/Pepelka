% Kako je bolje izracunati centre gruc?

N = 1000;
K = 3;

D = rand(N, 5);
labels = (randi(K,N,1));


[clusterIDs, ~, lbl] = unique(labels);
K = length(clusterIDs);

tic;
centers = zeros(K,size(D,2));
clsize = zeros(1,K);
for i = 1:K
    Dtmp = D(labels==i,:);
    centers(i,:) = mean(Dtmp,1);
    clsize(i) = size(Dtmp,1);
end
t1=toc

% compute binary membership matrix U
tic
E = eye(K);
U = logical(E(:,labels));
clsize2 = sum(U,2)';
centers2 = bsxfun(@rdivide, U*D, clsize');
t2 = toc


