% Test of DB index on Iris dataset 
%
% Ugotovitev: koda iz SOM toolboxa dela najbrž narobe, ker ne upošteva vseh
% možnih parov (gleda samo zgornji trikotnik v matriki R in ne celotnega).
% Èe bi upošteval vse možnosti (razen ko je i == j), bi dobil isto kot
% indexDB.

%% Load data
%load('iris_orig');


data = dlmread('iris.data','\t');
target = data(:,end);
data = data(:,1:end-1);


% Z-score of the data
%tmp = bsxfun(@minus,data,mean(data,1));
%data = bsxfun(@rdivide,tmp,std(data));


%% Run clusterer
% Run KM 10 times, pick result with best (max) SF score
Kmax = 12;
iters = 1;

val = zeros(1,Kmax-1);
val_som = zeros(1,Kmax-1);

for k=2:Kmax
    
    val_i = zeros(1,iters);
    val_i_som = zeros(1,iters);
    
    for r=1:iters
        %[labels,centroids,sumd] = kmeans(data,k, 'emptyaction', 'singleton');
        labels = pplk_runClusterer('AL',data,k,1);
        val_i(r) = indexDB(data,labels,0);
        val_i_som(r) = db_index(data, labels);
    end
    
    [val(k-1), bestInd ]= min(val_i); 
    [val_som(k-1), bestInd ]= min(val_i_som); 
end


%% Plot index values
% compute optimum
[OPT_val, OPT_k ]= min(val);
[OPT_val_som, OPT_k_som ]= min(val_som);


figure();
subplot(2,1,1); plot(2:Kmax, val,'.-'); hold on; plot(OPT_k+1,OPT_val,'ro'); title('DB code');
subplot(2,1,2); plot(2:Kmax, val_som,'.-'); hold on; plot(OPT_k_som+1,OPT_val_som,'ro'); title('SOM');


