 
load V_D31.mat

nClusts = 500;
Vd = V(:,1:nClusts);
N = size(Vd,1);
%%%%%%% Normalize rows of V
normV = sqrt(sum(Vd.^2,2));
Vd = bsxfun(@rdivide,Vd,normV+1e-10);

ticID = tic();
%%% Initialize centers
mu = zeros(nClusts,size(Vd,2));
% First center is set to one entry picked randomly
%[~,pp] = max(rand(size(Vd,1),1));
pp = 22;
mu(1,:) = Vd(pp,:);
% The other centers are selected to be farthest from previous centers
for i=2:nClusts
    ip = Vd*mu';
    minip = max(abs(ip'));
    [~,ii] = min(minip);
    mu(i,:) = Vd(ii,:);
end

toc(ticID)