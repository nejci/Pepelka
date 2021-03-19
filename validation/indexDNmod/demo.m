% Demonstration of the Modified Dunn's cluster validity index.
clc;
fprintf(1,'Cluster validation using modified Dunn''s index.\n');
fprintf(1,'--------------------------------------------------------\n');

%% Load dataset.
dataName = 'wave';
data = pplk_loadData(dataName);
[n,d]=size(data);
fprintf(1,'Dataset %s loaded!\nSamples: %d, features: %d.\n',dataName,n,d);
fprintf(1,'--------------------------------------------------------\n');

%% Cluster the data. 
% Utilize hierarchical single-linkage algorithm with Euclidean distance. 
% Number of clusters: from 2 to maxK.
% maxK = ceil(sqrt(n));
maxK = 9;
labels = zeros(n,maxK-1);

% Single-linkage clustering algorithm (from the Statistics toolbox).
fprintf(1,['Performing single-linkage clustering for k=2-',num2str(maxK),...
                                                          '  [']);
D=pdist(data,'euclidean');
htree=linkage(D,'single');

for k = 2:maxK
    % Store labels.
    labels(:,k-1)=cluster(htree,k);
    fprintf(1,'.');
end
fprintf(1,']\n');


% Plot clustering results.
fprintf(1,'Clustering results are shown in figure.\n');
fprintf(1,'--------------------------------------------------------\n');
options_SP.title=['Clusterings of ', dataName, ' with SL'];
options_SP.subtitle=num2cell(2:maxK);
options_SP.markerSize = 5;
pplk_scatterPlot(data,labels,[],options_SP);

%% Validate the clusterings with internal indices.
% Compare Modified Dunn's index to generalized version (Pal & Biswas, 1997)

% Construct a neighborhood graph on data points
% Choose type of graph: 'rng', 'gabriel', 'directedKnn', 'mutualKnn',
% 'symKnn', 'EMST' or 'epsilon'.
graph_type = 'gabriel';
isDirected = 0; % set to 1 when using directed graph, e.g. directedKnn
options.show = 0; % set it to 1 to enable plotting of a graph

fprintf(1,'Creating %s graph on data points ... ', graph_type);

tic();
[G,d,uniqueInd]=graph_create(data,[],graph_type,options);
toc()

if ~isempty(uniqueInd)
    data = data(uniqueInd,:);
    labels = labels(uniqueInd,:);
end
fprintf('done!\n');

fprintf(1,'Validating results with internal indices        [');

DNg = zeros(1,maxK-1);
DNs = zeros(1,maxK-1);

for k = 1:maxK-1
    %   -> generalized Dunn's index (Pal & Biswas, 1997)
    DNg(k) = indexDNg_graph(G,labels(:,k),data);
    
    %   -> modified Dunn's index (Ilc, 2012)
    DNs(k) = indexDNs_graph(G,data,labels(:,k),isDirected);
    
    % SLOWER: alternative calls, when graph G is not pre-computed
    %DNs(k) = indexDNg(data,labels(:,k),graph_type,options);
    %DNs(k) = indexDNs(data,labels(:,k),graph_type,options);
    
    fprintf(1,'.');
end
fprintf(1,']\n');

% find maximum index value
[val_g,ind_g]=max(DNg);
[val_s,ind_s]=max(DNs);

figure();
subplot(2,1,1); plot(2:maxK, DNg,'.-k');hold on; plot(ind_g+1,val_g,'sr');
title('Generalized Dunn''s index');
xlabel('number of clusters'); ylabel('index value');

subplot(2,1,2); plot(2:maxK, DNs,'.-k');hold on; plot(ind_s+1,val_s,'sr');
title('Modified Dunn''s index');
xlabel('number of clusters'); ylabel('index value');

fprintf(1,'Values of validation indices are shown in figure.\n');
fprintf(1,'Optimal number of clusters is indicated with red square.\n');
fprintf(1,'--------------------------------------------------------\n');

