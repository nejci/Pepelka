function [labels, moreInfo]=pplk_runClustererDel(methodName,data,K,nRuns,params)
% PPLK_RUNCLUSTERERDEL
% Run single clustering algorithm with data columns deletion.
% [labels,moreInfo] = pplk_runClustererDel(methodName,data,K,nRuns,params)
% -------------------------------------------------------------------------
% INPUTS
%   methodName  (string)  abbreviated name of clustering algorithm
%   data        (matrix)  input data, size of [N x D]
%               (string)  filename or name of data in the datasets folder
%   K           (int)     desired number of clusters
%   nRuns       (int)     number of runs of method
%   params      (struct)  optional structure with fields containing required
%                         parameter's value. Leave out for defaults;
%
% OUTPUTS
%   labels      (vector)  vector of data labels - clustering
%   moreInfo    (cell)    other info from clusterer (time, etc.)
%
%
% DESCRIPTION
%   Cluster data into K clusters with column deletion. Repeat clustering
%   D-times (D is number of features) and on each step skip one column of
%   data (on i-th step delete i-th column of data).
%
%   List of available methods with associated parameters:
%   --------------+-----------------------------------------+--------------
%   methodName    | params.methodName_                      | Description
%   --------------+-----------------------------------------+--------------
%   AL            | distance*                               | Average-linkage
%   CLUSOT        | cluster_method,msize**, shape, res,     | Clusot: SOM+RF/GB
%                 | theta, theta0 / g                       |
%   CS            | sigma**, Kin, Nin                       | Cauchy-Schwarz divergence
%   EM            | /                                       | Expectation Maximization
%   FCM           | /                                       | Fuzzy C-means
%   GSOM          | G,dG,alfa,maxIter, pOut, msize**, shape | Gravitational SOM
%                 | distance, showSOM, showGrav, advanced   |
%   HCL           | clustMethod, distance                   | Hierarchical Clustering (general, including AL,SL,WL and few others)
%   SOMKM         | nRuns, msize**, shape                   | K-means on the SOM
%   KCC           | /                                       | K-centers
%   KM            | maxIter, nRuns, distance                | K-means
%   KVV           | sigma*                                  | Kannan, Vempala and Vetta spectral
%   NC            | /                                       | Normalized cuts
%   NJW           | sigma*                                  | Ng, Jordan and Weiss spectral
%   RANDOM        | /                                       | random partition
%   SL            | distance*                               | Single-linkage
%   SPECLS        | Knn, mode                               | Spectral local scaling
%   WL            | distance*                               | Ward-linkage
%  ---------------+-----------------------------------------+--------------
%   *  can not-exists
%   ** can be []
%
% -------------------------------------------------------------------------
% Last modification: 27. September 2013
% (C) Pepelka Package, Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% -------------------------------------------------------------------------

callDir=chdir(pplk_homeDir());


%defaults
if ~exist('nRuns','var') || isempty(nRuns)
   nRuns = 1;
end
if ~exist('params','var') || isempty(params)
   params=pplk_setParamsDefault();
end


% CHANGE DIR - navigate to ..\methods\
%-------------------------------------------------------------------------
oldPath = chdir(['..',filesep,'methods']);



% USER DEFINED - read config file with user-defined functions
%-------------------------------------------------------------------------
f=fopen('methods.info','r');
usrfile=textscan(f,'%s %s %s %*[^\n]','delimiter','\t','commentStyle','%');
fclose(f);

usrMeth=usrfile{1};
usrAbbr=usrfile{2};
usrDesc=usrfile{3};



% RUN METHOD - wraper functions
%-------------------------------------------------------------------------
[N,D]=size(data);

labels=zeros(N,D,nRuns);
moreInfo=cell(1,nRuns); %to save additional information from method


fullMethodName=['pplk_clusterer',methodName];

% BUILT-IN functions -	filename starts with pplk_clusterer
%						have to return 2 output arguments: labels, moreInfo
if exist(fullMethodName,'file')
	for ind=1:nRuns
        for del=1:D
            cols = true(1,D);
            cols(del)=0; % skip del-th column
            [labels(:,del,ind), moreInfo{ind}]=feval(str2func(fullMethodName),data(:,cols),K,params);
        end
	end
	
% USER-DEFINED functions
else
	% Check for user-defined imported functions checking file
 	% ..\methods\methods.info .

	% search is made in the list of abbreviations, stored in
	% usrAbbr string array.
	idx = strmatch(upper(methodName),upper(usrAbbr),'exact');

	if ~isempty(idx) && exist(usrMeth{idx},'file')
		disp(['pplk_runClusterer :: Executing USER-DEFINED function ',usrAbbr{idx},': ',usrMeth{idx},'(...) -- ',usrDesc{idx}])
		
		% if user-defined method returns only one argument (labels)
		if nargout(usrMeth{idx}) == 1
            for ind=1:nRuns
                for del=1:D
                    cols = ones(1,D);
                    cols(del)=0; % skip del-th column
                    [labels(:,del,ind)]=feval(str2func(usrMeth{idx}),data(:,cols),K,params);
                end
            end
		else
            for ind=1:nRuns
                for del=1:D
                    cols = ones(1,D);
                    cols(del)=0; % skip del-th column
                    [labels(:,del,ind), moreInfo{ind}]=feval(str2func(usrMeth{idx}),data(:,cols),K,params);
                end
            end
		end
	else
		disp(['pplk_runClusterer :: Non-existing method name: ', methodName, '! Ignoring.']);
		labels=NaN;
	end
end

% CHANGE DIR - navigate back to the CORE folder
%-------------------------------------------------------------------------
chdir(oldPath);

%TODO
if params.showResult
    nRuns_show=nRuns;
	
	 % show first five of them
	if nRuns>5
        nRuns_show=5;
        disp('Only the first five results are displayed.');
	end
	
    for ind=1:nRuns_show
		if strcmp(methodName,'GSOM')
			% show auto-determined number of clusters
			K=moreInfo{ind}.nClust;
		end
        options.title=[methodName,', tek: ',num2str(ind)];
        options.colorMode='color';
        options.normalize=0;
        pplk_scatterPlot(data,labels(:,ind),K,options);
    end
end

chdir(callDir);