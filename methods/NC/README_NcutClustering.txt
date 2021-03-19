%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	       Normalized Cut Clustering Code                       %
%						                    %	
%  Timothee Cour (UPENN), Stella Yu (Berkeley), Jianbo Shi (UPENN)  %
%     						                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Installation Notes :

1) After you unzipped the files to mydir, 
   put the Current Directory in Matlab to mydir

2) In the matlab command prompt,

	type demoNcutClustering to see a demo

		or...
 
 	type main to initialize the paths to subfolders

3) You can now try any of the functions

The files were tested under matlab 6.5

Top level functions:

ncutW:  Given a similarity graph "W", computes Ncut clustering on the graph into "ncCluster" groups;
        NcutDiscrete = ncutW(W,nbCluster);

Use demo.m to see how to use the clustering function.
