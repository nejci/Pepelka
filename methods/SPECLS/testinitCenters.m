 
load V_D31.mat

nClusts = 500;
Vd = V(:,1:nClusts);
N = size(Vd,1);
%%%%%%% Normalize rows of V
normV = sqrt(sum(Vd.^2,2));
Vd = bsxfun(@rdivide,Vd,normV+1e-10);

ticID = tic();

%%%%%%%%%%%%%%%%%%%%%%%%%%
D = squareform(pdist(Vd,'euclidean'));
mu2 = zeros(nClusts,size(Vd,2));
Dtmp = zeros(nClusts,N);
p1 = 22; %randi(N,1);
mu2(1,:) = Vd(p1,:);
Dtmp(1,:) = D(p1,:);
for i=2:nClusts
   S = sum(Dtmp,1);
   [v,p2] = max(S);
   Dtmp(i,:) = D(p2,:);
   mu2(i,:) = Vd(p2,:);
end

toc(ticID)


