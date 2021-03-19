% function ClusterBasicsTest
%
% DESCRIPTION
%   conducts a quick test of the functions in the
%   CLUSTERBASICS package
% REFERENCE
%   please refer to the following paper if you use CLUSTERBASICS
%     A. Strehl, J. Ghosh and R. Mooney, "Impact of Similarity
%     Measures on Web-page Clustering", Proc. of the 17th National
%     Conference on Artificial Intelligence: Workshop of Artificial
%     Intelligence for Web Search (AAAI 2000), July 2000, Austin,
%     Texas, pp. 58-64
% RELEASE
%   version 1.0, 2002/04/20, tested on Matlab 5.2.0.3084, LNX86
%   available from http://www.strehl.com
%   license granted for research use ONLY (see README)
%   copyright (c) 1998-2002 by Alexander Strehl

function ClusterBasicsTest

% define some data matrix
x = [6.1 2.9 4.7 1.4;6.3 3.4 5.6 2.4;6 3 4.8 1.8;5.8 2.7 5.1 1.9;6.5 3.2 5.1 2;6.9 3.1 5.4 2.1;5 2 3.5 1;5 3.6 1.4 0.2;5.2 3.5 1.5 0.2;5.7 4.4 1.5 0.4;5.5 2.3 4 1.3;6.5 3 5.8 2.2;4.8 3 1.4 0.1;5.2 4.1 1.5 0.1;4.7 3.2 1.3 0.2;6.3 2.9 5.6 1.8;5.5 4.2 1.4 0.2;4.9 3.1 1.5 0.1;6.3 2.5 5 1.9;7.2 3 5.8 1.6;4.8 3.4 1.6 0.2;5.8 2.7 4.1 1;5 2.3 3.3 1;6 2.7 5.1 1.6;5.4 3.7 1.5 0.2;5.4 3.9 1.3 0.4;5.6 2.9 3.6 1.3;7.3 2.9 6.3 1.8;5.7 2.8 4.1 1.3;5.8 2.7 5.1 1.9;6.5 3 5.2 2;5.6 3 4.5 1.5;6.4 2.8 5.6 2.1;6.4 3.2 4.5 1.5;7.1 3 5.9 2.1;7.7 2.6 6.9 2.3;4.8 3 1.4 0.3;6.2 2.2 4.5 1.5;5.5 2.4 3.8 1.1;6.3 2.8 5.1 1.5;5.7 2.9 4.2 1.3;7.6 3 6.6 2.1;5 3.4 1.6 0.4;5.5 2.5 4 1.3;4.9 3 1.4 0.2;6.3 2.3 4.4 1.3;6.7 3 5.2 2.3;5.6 2.8 4.9 2;7.9 3.8 6.4 2;7.7 3 6.1 2.3;4.9 3.1 1.5 0.1;5.3 3.7 1.5 0.2;4.3 3 1.1 0.1;4.8 3.4 1.9 0.2;6.9 3.1 5.1 2.3;5.1 3.8 1.6 0.2;5.1 3.3 1.7 0.5;5.5 2.6 4.4 1.2;5.8 2.6 4 1.2;6.4 2.7 5.3 1.9;6.3 2.5 4.9 1.5;6.8 2.8 4.8 1.4;7.2 3.2 6 1.8;5.7 2.5 5 2;6.4 2.8 5.6 2.2;7 3.2 4.7 1.4;6.8 3.2 5.9 2.3;5.7 2.6 3.5 1;6.7 3.1 4.7 1.5;4.9 2.4 3.3 1;6.7 3.3 5.7 2.5;6.9 3.2 5.7 2.3;5.7 3.8 1.7 0.3;7.4 2.8 6.1 1.9;6 2.9 4.5 1.5;6.1 2.6 5.6 1.4;6.7 3.1 5.6 2.4;5.1 3.8 1.9 0.4;5.4 3 4.5 1.5;5.1 3.5 1.4 0.2;4.9 3.1 1.5 0.1;5 3.5 1.3 0.3;5.1 3.5 1.4 0.3;6.3 3.3 4.7 1.6;4.4 2.9 1.4 0.2;6.1 3 4.6 1.4;5.6 3 4.1 1.3;6.2 2.9 4.3 1.3;5 3.5 1.6 0.6;4.7 3.2 1.6 0.2;5.6 2.5 3.9 1.1;6.9 3.1 4.9 1.5;5.8 4 1.2 0.2;5.4 3.9 1.7 0.4;6 3.4 4.5 1.6;4.6 3.2 1.4 0.2;5.5 2.4 3.7 1;5.2 2.7 3.9 1.4;7.7 2.8 6.7 2;6.3 3.3 6 2.5;6.1 2.8 4 1.3;6.7 2.5 5.8 1.8;6.4 2.9 4.3 1.3;5.7 2.8 4.5 1.3;6.1 2.8 4.7 1.2;5.4 3.4 1.7 0.2;5.8 2.8 5.1 2.4;4.9 2.5 4.5 1.7;4.5 2.3 1.3 0.3;6.6 3 4.4 1.4;6.8 3 5.5 2.1;4.4 3 1.3 0.2;6 2.2 4 1;4.6 3.4 1.4 0.3;5.2 3.4 1.4 0.2;6.4 3.2 5.3 2.3;6.7 3.1 4.4 1.4;5.8 2.7 3.9 1.2;5.9 3 4.2 1.5;5.5 3.5 1.3 0.2;6.4 3.1 5.5 1.8;6.2 3.4 5.4 2.3;6.5 2.8 4.6 1.5;6.6 2.9 4.6 1.3;5.6 2.7 4.2 1.3;5.1 3.8 1.5 0.3;5 3.4 1.5 0.2;7.2 3.6 6.1 2.5;5.7 3 4.2 1.2;4.6 3.6 1 0.2;5.4 3.4 1.5 0.4;5.9 3 5.1 1.8;5.1 3.7 1.5 0.4;6.2 2.8 4.8 1.8;5 3.3 1.4 0.2;5.1 3.4 1.5 0.2;5.9 3.2 4.8 1.8;6.3 2.7 4.9 1.8;4.8 3.1 1.6 0.2;5 3 1.6 0.2;4.6 3.1 1.5 0.2;5.1 2.5 3 1.1;6.7 3.3 5.7 2.1;7.7 3.8 6.7 2.2;4.4 3.2 1.3 0.2;6.1 3 4.9 1.8;6.5 3 5.5 1.8;5 3.2 1.2 0.2;6 2.2 5 1.5;6.7 3 5 1.7];

% and some original labeling
cat = [2 3 3 3 3 3 2 1 1 1 2 3 1 1 1 3 1 1 3 3 1 2 2 2 1 1 2 3 2 3 3 2 3 2 3 3 1 2 2 3 2 3 1 2 1 2 3 3 3 3 1 1 1 1 3 1 1 2 2 3 2 2 3 3 3 2 3 2 2 2 3 3 1 3 2 3 3 1 2 1 1 1 1 2 1 2 2 2 1 1 2 2 1 1 2 1 2 2 3 3 2 3 2 2 2 1 3 3 1 2 3 1 2 1 1 3 2 2 2 1 3 3 2 2 2 1 1 3 2 1 1 3 1 3 1 1 2 3 1 1 1 2 3 3 1 3 3 1 3 2];

draw_clustering2(x,3,cat,'Original',1,0)

pause

% conduct a variety of clusterings using different algorithms and similarity measures
sim='simeucl';

%%
cl1 = clgraph(x,3,sim); %weighted graph partitioning using METIS
draw_clustering2(x,3,cl1,'1. Graph',1,0);


% evaluate the clusterings
BAL=evalbalance([],cl1,[],[]); % balance of clustering
MSE=evalmse([],cl1,x,sim); %Mean-squared error
F1=evalf(cat,cl1,[],[]); %f1 measure
NMI=evalmutual(cat,cl1,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

%%
cl2 = clkmeans(x,3,sim); % kmeans
draw_clustering2(x,3,cl2,'2. Kmeans',1,0);

% evaluate the clusterings
BAL=evalbalance([],cl2,[],[]); % balance of clustering
MSE=evalmse([],cl2,x,sim); %Mean-squared error
F1=evalf(cat,cl2,[],[]); %f1 measure
NMI=evalmutual(cat,cl2,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

%%
cl3 = clcgraph(x,3,sim); % weighted balanced graph partitioning using PMETIS
draw_clustering2(x,4,cl3,'3. Cgraph- pmetis',1,0);

% evaluate the clusterings
BAL=evalbalance([],cl3,[],[]); % balance of clustering
MSE=evalmse([],cl3,x,sim); %Mean-squared error
F1=evalf(cat,cl3,[],[]); %f1 measure
NMI=evalmutual(cat,cl3,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

%%
cl4 = clagglmin(x,3,sim); % single-linkage (slow implementation)
draw_clustering2(x,3,cl4,'4. Single-linkage',1,0);

% evaluate the clusterings
BAL=evalbalance([],cl4,[],[]); % balance of clustering
MSE=evalmse([],cl4,x,sim); %Mean-squared error
F1=evalf(cat,cl4,[],[]); %f1 measure
NMI=evalmutual(cat,cl4,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

%%
cl5 = clrand(x,3,[]); % random clustering
draw_clustering2(x,3,cl5,'5. Random',1,0);

% evaluate the clusterings
BAL=evalbalance([],cl5,[],[]); % balance of clustering
MSE=evalmse([],cl5,x,sim); %Mean-squared error
F1=evalf(cat,cl5,[],[]); %f1 measure
NMI=evalmutual(cat,cl5,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

%%
cl6 = clhgraph(x,3,[]); % hypergraph partitioning using HMETIS
draw_clustering2(x,3,cl6,'6. Hmetis',1,0);

% evaluate the clusterings
BAL=evalbalance([],cl6,[],[]); % balance of clustering
MSE=evalmse([],cl6,x,sim); %Mean-squared error
F1=evalf(cat,cl6,[],[]); %f1 measure
NMI=evalmutual(cat,cl6,[],[]); %NMI measure

fprintf('BAL: %f \nMSE: %f\nF1: %f\nNMI: %f\n', BAL, MSE, F1, NMI);
%%
pause

close all;

