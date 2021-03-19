clear

% load the data and initial model parameters
%
%   k = number of ensemble clusters
%   N = number of base clusterings
%   M = number of data points
%
% base_labels:              M*N, base clustering results to be processed
% Palpha:                   k*1, initial value for the model parameter of Dirichlet distribution
% Pbeta:                    cell with N elements for N base clusterings, each element is a
%                           k*q matrix, i.e., initial value for k parameters of a q-dimensional discrete distribution
% number_baseclusterers:    1*N, number of clusters in each of N base clustering results
load Iris.mat;

% PramaLap:                 parameter for laplace smoothing
PramaLap=0.000001 ;

% If use random initialization
Kcons = length(unique(true_labels));
N = size(base_labels,2);

Palpha = rand(Kcons,1);
Pbeta = cell(1,N);
for i=1:length(Pbeta)
    q = length(unique(base_labels(:,i)));
    temp = rand(Kcons,q);
    temp = temp./(sum(temp,2)*ones(1,q));    
    Pbeta{i} = temp;
end



% learn BCE 
[phiAll,gammaAll,resultAlpha,resultBeta]=learnBCE(base_labels,Palpha,Pbeta,PramaLap,number_baseclusterers);


% calculate accuracy
k=length(unique(true_labels));
M=length(true_labels);

% Obtain the cluster assignments from BCE
Ensemble_labels=zeros(1,M);
for index=1:M
    wtheta(:,index)=gammaAll(:,index);
    bb=find(wtheta(:,index)==max(wtheta(:,index)));
    Ensemble_labels(index)=bb(1);
end

% Calculate the accuracy based on true labels and BCE results
[accu]=calculateAccuracy(true_labels,Ensemble_labels);
disp(['The micro-precision of BCE is ',num2str(accu)]);




