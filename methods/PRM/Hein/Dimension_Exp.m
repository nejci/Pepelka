function Dimension_Exp(set)
% This function reproduces the experimental results from the paper
% M. Hein, J.-Y. Audibert, Intrinsic dimensionality estimation of
% Submanifolds in R^d, Proceedings of the 22nd ICML, 289-296, Eds. L. de
% Raedt and S. Wrobel, 2005
%
% Usage: Dimension_Exp(set)
% where the possible sets are
% 1: Sinusoid
% 2: Sphere
% 3: Gaussian distribution with full support
% 4: Moebius strip
% 5: 12-dimensional submanifold of R^72

path2data = ['data', filesep];
datasets={'Sinusoid','S','Gauss','Moebius','M12'};
dim_right={[1],[3,5,7,9],[3,4,5,6],[2],[12]};
dim={[3],[4,6,8,10],[3,4,5,6],[3],[72]};
num={[400,500,600], [600, 800, 1000, 1200], [100, 200, 400, 800], [20, 40, 80, 120], [200, 400, 800, 1600]};
runs=90;
est=zeros(size(dim{1,set},2),size(num{1,set},2),runs ,3);

for k=1:size(dim{1,set},2)
  for i=1:size(num{1,set},2)
    if(size(dim{1,set},2)>1)
      name=[path2data,datasets{set},filesep,datasets{set},num2str(dim_right{1,set}(k)),'_',num2str(num{1,set}(i)),'.BIN'];
    else
      name=[path2data,datasets{set},filesep,datasets{set},'_',num2str(num{1,set}(i)),'.BIN'];  
    end
    fid=fopen(name,'rb');
    X=fread(fid,inf,'float32');
    fclose(fid);
    X=reshape(X,dim{1,set}(k),runs*num{1,set}(i));
    for j=1:runs
       Z=X(:,(j-1)*num{1,set}(i)+1:j*num{1,set}(i));
       est(k,i,j,:)=GetDim(Z);
    end
  end
end

for k=1:size(dim{1,set},2)
  for i=1:size(num{1,set},2)
    disp(['Dimension of the ambient space: ',num2str(dim{1,set}(k)),' Correct intrinsic dimension: ',num2str(dim_right{1,set}(k)),', Number of Points: ',num2str(num{1,set}(i))])
    disp(['Correct estimates of the Int. Dim. Estimator  : ', num2str(sum(est(k,i,:,1)==dim_right{1,set}(k))),' out of ',num2str(runs)])
    disp(['Correct estimates of the Correlation Dimension: ', num2str(sum(abs(est(k,i,:,2)-dim_right{1,set}(k))<0.5)),' out of ',num2str(runs)])
    disp(['Correct estimates of the Takens Estimator     : ', num2str(sum(abs(est(k,i,:,3)-dim_right{1,set}(k))<0.5)),' out of ',num2str(runs)])
  end
end
    