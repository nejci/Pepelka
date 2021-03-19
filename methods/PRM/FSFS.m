function [redu,fwt]=FSFS(data,k,method)
%FSFSPACK Reduces a feature set with Feature Selection using Feature Similarity 
% [redu,fwt]=FSFS(data,original_size,k)
%
% INPUTS
% data   : Data matrix containing the original feature set. Each column
%          represents a feature.
% k      : Scale parameter which decides the size of the reduced feature set.
%          Approximately, k = #original features - #reduced feature set.
% method : method for calculating similarity between features:
%           - 'cor' = Correlation Coefficient
%           - 'lin' = Linear Regression error
%           - 'mic' = Maximal Information Compression Index (default)
%
% OUTPUTS
% redu : Reduced features. A vector containing the feature numbers of the 
%        original feature set which are retained in the reduced set.
% fwt :  feature weights of the features in redu. 
%
% Source written by Pabitra Mitra.
% Code clean-up and function interface adaptation made by Nejc Ilc.
%
% Reference:
% P. Mitra, C. A. Murthy and S. K. Pal, Unsupervised Feature Selection using 
% Feature Similarity, IEEE Transactions on Pattern Analysis and Machine
% Intelligence, Vol. 24, No. 4, pp 301-312, April 2002.

[no_data, no_feature]=size(data);
kk=[];

if ~exist('k','var') || isempty(k)
    error('Parameter k is not defined!');
end 

if ~exist('method','var') || isempty(method)
   warning('pplk:FSFS','No method specified, using "Maximal Information Compression Index" as default.');
   method = 'mic';
end

method = lower(method);


if k >= no_feature
    error('k should be less than number of features.');
end



switch method
    % 1 = Feature Similarity: Correlation Coeff
    % 2 = Feature Similarity: Linear Regression error
    % 3 = Feature Similarity: Maximal Information Compression Index
    case 'cor'
        method = 1;
    case 'lin'
        method = 2;
    case 'mic'
        method = 3;
    otherwise
        error('Wrong method specified!');
end

% Form the inter-feature distance matrix (Upper Triangular)
%fprintf(1,'Computing Feature Similarities..\n');
dm = zeros(no_feature);

for i = 1 : (no_feature-1)
   %fprintf(1,'Similarity Computed for Feature No. %d\n',i);
   for j = (i+1) : no_feature
      x1=data(:,i);
      x2=data(:,j);
      dm(i,j)=f2f(x1,x2,method);
   end
end

drift=1.0;
% Form a vector containing the distance of k-NN for each feature.
kd = zeros(1,no_feature);
for i=1:no_feature,
   if i==1
      dd=dm(1,2:no_feature);
   elseif i==no_feature
      dd=dm(1:no_feature-1,no_feature)';
   else
      dd=[dm(i,i+1:no_feature),dm(1:i-1,i)'];
   end
   dd=sort(dd);
   kd(i)=dd(k);
end
kd0=kd; % Store the original r_k's


% Condense the feature set
%fprintf(1,'\nClustering the Features..\n');
rfs=[];
rfd=[];
ee=[];
low=9999;
iter=0;
prev_lower=9999;
tagm=ones(1,no_feature);
while (no_feature-trace(dm)) >0,
   iter=iter+1;
   if k > (no_feature-trace(dm)-1)
      k=no_feature-trace(dm)-1;
   end
   if k<=0
      break;
   end
   prev_lower=low;
   [low,fetr]=lowerror(dm,k);
  
   
   % Adjust k
   while low > drift*prev_lower,
      k=k-1;
      if k==0
         break;
      end
      [low,fetr]=lowerror(dm,k);
   end
   
   if k <=0
      break;
   end
   dm=updatedm(dm,fetr,k);
          
   kk=[kk;k];
   ee=[ee;low];
     
   tagm(fetr)=0;
   for i=1:no_feature,
      for j=1:no_feature,
         if dm(i,i)==1
            tagm(i)=0;
         end
      end
   end
   
end

for i=1:no_feature,
   if dm(i,i)==0
      rfs=[rfs;i];
      rfd=[rfd;kd0(i)];
   end
end

%fprintf(1,'Features Clustered.\n');
     
redu=rfs;
fwt=rfd;

end

% Helper functions

function dist=f2f(x1,x2,method)
% Function to compute correlation between two variables

no_x1=size(x1,1);
no_x2=size(x2,1);
dist=0.0;

% Distance = Correlation Cefficient
if method==1
    num=0.0;den1=0.0;den2=0.0;
    x1bar=mean(x1);x2bar=mean(x2);
    for i=1:no_x1,
        num=num+abs(x1(i)*x2(i)-x1bar*x2bar);
        den1=den1+(x1(i)-x1bar)^2;den2=den2+(x2(i)-x2bar)^2;
    end
    dist=num/sqrt(den1*den2);
    
    
% Distance = linear regression error
elseif method==2
    num=0.0;den=0.0;x1bar=mean(x1);x2bar=mean(x2);
    for i=1:no_x1,
        num=num+x1(i)*x2(i)-x1bar*x2bar;
        den=den+x1(i)^2-x1bar^2;
    end
    a=num/den;b=x2bar-a*x1bar;
    
    dist=0.0;
    
    for i=1:no_x1,
        dist=dist+abs(x2(i)-a*x1(i)-b)/sqrt(a^2+b^2);
    end
    dist=dist/no_x1;

% Distance = Maximum information compression index
elseif method==3
    sxy=0.0;sx=0.0;sy=0.0;mnx1=0.0;mnx2=0.0;
    for i=1:no_x1,
        sxy=sxy+x1(i)*x2(i);
        sx=sx+x1(i)^2;
        sy=sy+x2(i)^2;
        mnx1=mnx1+x1(i);
        mnx2=mnx2+x2(i);
    end
    mnx1=mnx1/no_x1;
    mnx2=mnx2/no_x2;
    
    sxy=(sxy/no_x1)- mnx1*mnx2;
    sx=(sx/no_x1)-mnx1^2;
    sy=(sy/no_x1)-mnx2^2;
    
    if (sx-sy) ==0
        theta=0.5*pi/2;
    else
        theta=0.5*atan(2*sxy/(sx-sy));
    end
    
    a=-cot(theta);
    b=mean(x1)*cot(theta)+mean(x2);
    
    dist=0.0;
    for i=1:no_x1,
        dist=dist+abs(x2(i)-a*x1(i)-b)/sqrt(a^2+b^2);
    end
    dist=dist/no_x1;
    
elseif method==4
    cv=cov(x1,x2)/no_x1;
    dist=min(eig(cv));
end

end

function dm1=updatedm(dm,indx,k)
% Function to recompute the distance matrix during clustering

no_feature=size(dm,1);
HIGH=9999;
i=indx;

if i==1
    dd=[HIGH,dm(1,2:no_feature)];
    for l=1:no_feature,
        if dm(l,l)==1
            dd(l)=HIGH;
        end
    end
    
elseif i==no_feature
    dd=[dm(1:no_feature-1,no_feature)',HIGH];
    for l=1:no_feature,
        if dm(l,l)==1
            dd(l)=HIGH;
        end
    end
    
else
    dd=[dm(1:i-1,i)',HIGH,dm(i,i+1:no_feature)];
    for l=1:no_feature,
        if dm(l,l)==1
            dd(l)=HIGH;
        end
    end
    
end

[dd,dindx]=sort(dd);

dm1=dm;

for l=1:k,
    indx1=dindx(l);
    dm1(indx1,indx1)=1;
end
end

function [erro,indx]=lowerror(dm,k)
% Function to find the lowest error
no_feature=size(dm,1);
HIGH=9999;

for i=1:no_feature,
   if dm(i,i)==1
      kd(i)=HIGH;
   else
      
   if i==1
      dd=[HIGH,dm(1,2:no_feature)];
   elseif i==no_feature
      dd=[dm(1:no_feature-1,no_feature)',HIGH];
   else
      dd=[dm(1:i-1,i)',HIGH,dm(i,i+1:no_feature)];
   end
   for l=1:no_feature,
      if dm(l,l)==1
         dd(l)=HIGH;
      end
   end
   
   dd=sort(dd);
   kd(i)=dd(k);
   end

end

[erro,indx]=min(kd);
end

