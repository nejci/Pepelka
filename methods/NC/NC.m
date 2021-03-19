function [retVal]=NC(data,nbCluster,options)
% demoNcutClustering
% 
% demo for NcutClustering
% also initialize matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
% Modified by Nejc Ilc, 2013.


% deprecated - initialize matlab paths to subfolders - there are none now
% ADDED_PATH=main();

if ~exist('options','var')
    %parameters
    options.scaleSigma = [];
    options.offset = 0.5; %offset in the diagonal of W, default 0.5
    options.verbose = 1; %0 for verbose off mode, 1,2,3 for verbose on modes
    options.maxiterations = 100; %max number of iterations in eigensolver
    options.eigsErrorTolerance = 1e-8; %error tolerance in eigensolver
    options.valeurMin=1e-6; %truncates any values in W less than valeurMin
end


% make up a point data set
% caseid = 4;
% [data,size_cluster] = build_scene(caseid);

if options.verbose
    figure(1);clf;
    plot(data(1,:),data(2,:),'ks', 'MarkerFaceColor','k','MarkerSize',5); axis image; hold on; 
end
%disp('This is the input data points to be clustered, press Enter to continue...');
%pause;

if options.verbose
    disp('Compute clustering...');
end

ticID = tic();
% compute similarity matrix
W = compute_relation(data,options.scaleSigma);

% Cluster
NcutDiscrete = ncutW(W,nbCluster,options);
time = toc(ticID);

if options.verbose
	disp(['The computation took ' num2str(time) ' seconds']);
end

%transform n-outOf-1 presentation
labels=NcutDiscrete(:,1);
for labInd=2:nbCluster
    labels(NcutDiscrete(:,labInd)==1)=labInd;
end

if options.verbose
% display clustering result
    cluster_color = 'rgbmyc';
    figure(2);clf;
    for j=1:nbCluster,
        id = find(NcutDiscrete(:,j));
        plot(data(1,id),data(2,id),[cluster_color(j),'s'], 'MarkerFaceColor',cluster_color(j),'MarkerSize',5); hold on; 
    end
    hold off; axis image;
end

retVal.target=labels;
retVal.time=time;
retVal.options=options;


