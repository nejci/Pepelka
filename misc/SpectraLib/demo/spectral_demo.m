% Demonstrates the usage of some of the main functions in the spectral
% library. 

% read in the similarity and number of clusters 
[S,ca,points]=read_from_data_file('four-clouds'); 

% perform various clustering algorithms on it 
sl=single_linkage(S,4); 
mcut=mcut_kmeans(S,4);

% Show the results 

figure(1); 
plot2Dpoints_with_clusters(points,sl); 
title('Single Linkage on four-clouds'); 

figure(2); 
plot2Dpoints_with_clusters(points,mcut); 
title('Multicut-Kmeans  on four-clouds'); 

