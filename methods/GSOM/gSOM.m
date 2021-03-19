function [retVal]=gSOM(data,G,dG,alfa,maxIter,pOut,msize,shape,showSOM,showGrav,isMass,kPredef,distance,advanced)
%
%	E X P E R I M E N T A L    version!
%
% [retVal]=gSOM(data,G,dG,alfa,maxIter,pOut,msize,shape,showSOM,showGrav,isMass,kPredef,distance,advanced)
%--------------------------------------------------------------------------
% gSOM algorithm for clustering
% including some experimental features and a lot of parameters ;)
%
% Run gSOM_demo for an illustrative example on how to use the provided functions.
%--------------------------------------------------------------------------
% INPUTS
%   data:		nD x d matrix of data, where nD is the number of data points
%				and d the dimensionality of data points.
%	G:			gravitational coefficient [default=0.0008]
%	dG:			reduction factor of G at every iteration [deafult=0.045]
%	alfa:		merging distance
%	maxIter:	maximum number of iterations [default=100, 0=unlimited]
%	pOut:		probability of choosing random particle y instead of
%				strictly a neighbour of particle x. pOut is 1-p, where p is
%				the parameter described in [1]. [default=0.1]
%   msize:		size of map M=a x b (1x2 vector of a and b), e.g. [17,12];
%				if empty [], heuristic is used M=ceil(5*sqrt(nD)), where nD is the
%				number of data points;
%				map size can also be determined with provided scale factor S.
%				Toolbox recognizes this by negative sign, e.g. use msize=[-0.75] to
%				apply size of 0.75*5*sqrt(nD). Default is -1.
%	shape:		shape of the SOM grid ('rect' for rectangular - default, 'hexa' for
%				hexagonal)
%	showSOM:	if 1, trained SOM is displayed
%	showGrav:	if 1, gravitational clustering is visualized in detail
%				if 2, only the final result is plotted
%	isMass:		[experimental!]
%				if 1 or 2, mass of neurons is computed as number of data points that
%				belong to certain winning neuron. [default=0].
%				if 1, gravity force is proportional to the mass of neurons
%				if 2, gravity force is inverse proportional to the mass
%	kPredef:	optional parameter to force at least kPredef clusters in data [default=2]
%	distance:	which distance metric to use:
%					- 'sqEuclidean':  squared Euclidean distance (default)
%					- 'correlation':  one minus (centered) Pearson correlation between points
%					- 'ucorrelation': one minus (uncentered) Pearson correlation between points
%	advanced:	additional settings for training, given as name-value cell array.
%				Options marked with * are parameters for gravitational
%				algorithm. Other are SOM's parameters.
%               Options:
%					- 'algorithm'	('seq' | 'batch')
%					- 'radius_ini'	(numeric)
%					- 'radius_fin'	(numeric)
%					- 'alpha_ini'	(numeric)
%					- 'alpha_type'	('linear' | 'inv' | 'power')
%					- 'initalg'		('randinit' | 'lininit')
%					- 'trainlen'	('short' | 'long' | [numRoughPhase, numFinePhase])
%					- 'sampleOrder'	('ordered' | 'random')
%					- 'neigh'		('bubble' | 'gaussian' | 'cutgauss' | 'ep')
%                 * - 'saveItersState' (0 | 1) if 1, save state of iters
%				  * - 'eps'			('numeric') moving threshold - stopping criterion
%				  * - 'normalize'	(0 | 1) normalize data on [0 1] or on unit sphere
%				  * - 'dist_mode'	('euc_AR' | 'euc_pure' | 'euc_semi' |
%									 'dotprod_linAR' | 'dotprod_sqAR' |
%									 'dotprod_linA' | 'dotprod_sqA')
%					!! VALID ONLY WHEN distance == {'correlation','ucorrelation'} !!
%							'euc_AR':
%									force of gravity F and movement are computed
%									in Euclidean space. Additionaly, we compute
%									correlation r between two particles and use it
%									to weight the movement. r ranges from -1 to 1,
%									so the movement can be repulsive (-1..0) or
%									attractive (0..1). SOM is trained using
%									the same rule as for sqEuclidean distance.
%							'euc_semi':
%									correlation is used to compute distance
%									between two particles and euclidean
%									norm to compute the lenght of move
%									(norm(F)). SOM is trained using
%									the same rule as for sqEuclidean distance.
%							'euc_pure':
%									all distances are measured as
%									correlation. Drawback: unwanted behaviour may occur when
%									computing movement and there is a zero
%									move distance computed (d(x,x+F)) although
%									the force is positive. SOM is trained using
%									the same rule as for sqEuclidean distance.
%
%							'dotprod_linAR':
%									SOM is trained using dot-product rule.
%									Distance is measured as d=2*(1-r), where
%									r = x * y'. F=G*m1*m2/d.
%									If r > 0, F is attractive,
%									else if r < 0, F is repulsive.
%							'dotprod_linA':
%									SOM is trained using dot-product rule.
%									Distance is measured as d=2*(1-r), where
%									r = x * y'. F=G*m1*m2/d.
%									F is always attractive.
%							'dotprod_sqAR':
%									SOM is trained using dot-product rule.
%									Distance is measured as d=2*(1-r), where
%									r = x * y'. F=G*m1*m2/d^2.
%									If r > 0, F is attractive,
%									else if r < 0, F is repulsive.
%							'dotprod_sqA':
%									SOM is trained using dot-product rule.
%									Distance is measured as d=2*(1-r), where
%									r = x * y'. F=G*m1*m2/d^2.
%									F is always attractive.
%--------------------------------------------------------------------------
% OUTPUTS:
%   retVal :   return stucture with following fields
%				.target				: labels of data points - final clustering
%				.nClust				: number of final clusters
%				.iter				: number of iterations
%				.time				: elapsed time [SOM, Gravity, total]
%				.SOMprop			: properties of SOM
%				.SOMtrain			: info about SOM training
%				.kPredefReached		: flag that is 1 when gSOM stops
%									  because of reaching the predefined number of clusters
%
%
%------- LEGAL NOTICE -----------------------------------------------------
% gSOM - clustering method based on gravitational Self-Organizing Map
%
%	Copyright (C) 2012  Nejc Ilc
%
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	any later version.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	You should have received a copy of the GNU General Public License
%	(Licence-GPL.txt) along with this program.
%	If not, see: http://www.gnu.org/licenses/.
%
%------- REFERENCE --------------------------------------------------------
% N. Ilc, A. Dobnikar, Gravitational clustering of the self-organizing map,
% in: A. Dobnikar, U. Lotric, B. Ster (Eds.), Adaptive and Natural
% Computing Algorithms, Vol. 6594 of Lecture Notes in Computer Science,
% Springer Berlin / Heidelberg, 2011, pp. 11-20.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.2.1
% Last modified: 21-May-2014 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Bug reports, comments and questions are appreciated.
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================


%------- LOAD SOM TOOLBOX -------------------------------------------------
% Add a path to SOM Toolbox, modified to faster compute PCA, bug relating
% size determination of the SOM was smashed, some integration issues
% resolved (function som_make is now fully customizable, taking additional
% parameters through the cell array 'advanced').
addpath('..\..\misc\somtoolbox_pplk');

%------- CHECK-IN ---------------------------------------------------------

if ~exist('G','var') || isempty(G)
    G = 0.0008;
end
if ~exist('dG','var') || isempty(dG)
    dG = 0.045;
end
if ~exist('alfa','var') || isempty(alfa)
    alfa = 0.01;
end
if ~exist('maxIter','var') || isempty(maxIter)
    maxIter = 100;
end
if ~exist('pOut','var') || isempty(pOut)
    pOut = 0.1;
end
if ~exist('msize','var') || isempty(msize)
    msize = -1;
end
if ~exist('shape','var') || isempty(shape)
    shape = 'rect';
end
if ~exist('kPredef','var') || isempty(kPredef)
    kPredef=2;
end
if ~exist('isMass','var') || isempty(isMass)
    isMass=0;
end
if ~exist('showSOM','var') || isempty(showSOM)
    showSOM=0;
end
if ~exist('showGrav','var') || isempty(showGrav)
    showGrav=0;
end
if ~exist('distance','var') || isempty(distance)
    distance = 'sqEuclidean';
end

% shape can be also 1 for 'rect' or 2 for 'hexa'
if isnumeric(shape)
    if shape == 1
        shape = 'rect';
    elseif shape == 2
        shape = 'hexa';
    end
end

% advanced settings
% want to to store state of every iteration for postprocessing?
saveItersState = 0;
% threshold to determine when particles stop moving - a stopping criterion
epsStop = 1e-3;
normalize = 0;
%dist_mode = 'default';
som_mode = 'euc';
grav_mode = 'pure';
                        
if ~exist('advanced','var') || isempty(advanced)
    advanced={};
    
else
    % parse those related with gravitational algorithm
    i=1;
    while i<=length(advanced)
        if ischar(advanced{i})
            switch advanced{i}
                case 'saveItersState', i=i+1; saveItersState = advanced{i};
                case 'eps', i=i+1; epsStop = advanced{i};
                case 'normalize', i=i+1; normalize = advanced{i};
                case 'dist_mode', i=i+1; dist_mode=advanced{i};
                    if any(strcmp(dist_mode,...
                            {'euc_AR','euc_semi','euc_pure',...
                            'dotprod_linA','dotprod_sqA','dotprod_linAR','dotprod_sqAR'}))
                        [som_mode,grav_mode]=strtok(dist_mode,'_');
                        grav_mode=grav_mode(2:end);
                    else
                        warning(['Parameter dist_mode ''',dist_mode,''' is wrong! Using defaults.']);
                        dist_mode = 'default';
                    end
                    if strcmp(dist_mode,'default')
                        som_mode = 'euc';
                        grav_mode = 'pure';
                        dist_mode = 'euc_pure';
                        advanced{i}=dist_mode;
                    end
                    
                otherwise
                    i=i+1;
            end
        end
        i=i+1;
    end
end




%------- PREPROCESING -----------------------------------------------------


% Data normalization on interval [0,1], only when Squared Euclidean
% distance used - is it OK? Correlation implicitely normalizes data

dataStart=data;

if normalize
    switch distance
        
        % normalize - divide with the norm (to produce unit vectors)
        case 'ucorrelation'
            for ind=1:size(data,1)
                data(ind,:)= data(ind,:)./norm(data(ind,:));
            end
            
            % mean-shift (subtract mean) and normalize (divide with the norm)
        case 'correlation'
            for ind=1:size(data,1)
                datum=data(ind,:)-mean(data(ind,:));
                data(ind,:)= datum./norm(datum);
            end
            
            % rescale on [0,1] proportionaly
        case 'sqEuclidean'
            data = som_normalize(data, 'propor');
            
    end
end

[nPstart,dimStart]=size(dataStart);

%------- LEVEL ONE --------------------------------------------------------

ticID = tic();

% data representation with SOM neurons
[winners,sLabels,sMap] = createSOM(data,msize,shape,showSOM,isMass,distance,advanced);

tSOM=toc(ticID);

% from now on, data are considered as winning neurons
dataGrav = winners.coords;
len = size(dataGrav,1);

% ---------------------- danger --------------
% If we want to use parameters (G, alfa) in approximately the same order of
% magnitude for all data sets, we have to normalize data before applying
% gravitational algorithm. For Euclidean distance the most appropriate
% normalization is proportional scaling (to preserve the shape of dataset).
% What about when distance is correlation? Codebook vectors are already normalized!
if ~normalize && strcmp(distance,'sqEuclidean')
    %dataGrav = som_normalize(dataGrav, 'propor');
end

% analyze dataGrav
% DEBUG **********************************************
%DIST = pdist(dataGrav,'euclidean');
%fprintf(1,'gSOM debug | maxPdist: %f, minPdist: %f, meanPdist: %f, sumPdist: %f\n',max(DIST),min(DIST),mean(DIST),sum(DIST));
%fprintf(1,'SOM\t%f\t%f\t%f\t%f\n',max(DIST),min(DIST),mean(DIST),sum(DIST));
% normalize - bound max distance between points to 1
% find max distance (between points d1 and d2)
if strcmp(distance,'sqEuclidean')
    DIST = squareform(pdist(dataGrav,'euclidean'));
    maxDi_val = max(DIST,[],1);
    [maxD,maxD_indJ] = max(maxDi_val);
    d1 = maxD_indJ;
    % move all data, so that d1 is in the origin (0,0)
    dataNorm = bsxfun(@minus,dataGrav,dataGrav(d1,:));
    % divide all data points by maxD
    dataNorm = bsxfun(@rdivide, dataNorm, maxD);
    % move back to original position of d1
    dataGrav = bsxfun(@plus,dataNorm,dataGrav(d1,:));
end
% DEBUG **********************************************


% neighbourhood info
neighC = winners.neigh;

% mass of each neuron; if isMass==0, mass is a vector of ones
mass = winners.mass;

% showtime
if(showGrav>0)
    figure();
    if dimStart > 2
        %compute PCA
        if isempty(winners.dataPca)
            eig_vec=princomp(data,'econ');
            dataStartOK=(eig_vec(:,1:2)' * data')';
        else
            dataStartOK=winners.dataPca;
        end
        if isempty(winners.pca_vec)
            suns_pca_vec=eig_vec;
        else
            suns_pca_vec=winners.pca_vec;
        end
    else
        if normalize
            dataStartOK=data;
        else
            dataStartOK=dataStart;
        end
    end
    % if you want to see underlying data, uncomment next statement
    % plot(dataStartOK(:,1),dataStartOK(:,2),'bo');
    axis('equal')
end



%------- LEVEL TWO --------------------------------------------------------

% C is a cell of all particles that is updated on every merge of two
% particles
C = num2cell(1:len);

% list of indeces of remaining (valid) particles
valid = 1:len;

% curent number of clusters
nClust = len;

% flag to identify the case, when gSOM stopped because of predefined number
% of clusters
retVal.kPredefReached = 0;

% is the predefined number of clusters already reached?
if nClust <= kPredef
    retVal.kPredefReached=1;
end

if saveItersState == 1
    storedIters = [];
    storedIters.header = {'nClust', 'G', 'dataGrav','C','neighC','mass','valid'};
    storedIters.data = cell(1,7); 
end


tidID = tic();



% TODO - infinite loop

% MAIN LOOP over iterations
% Stopping criteria:
%   - maxIters reached, except if maxIter = 0
%   - movements of particles are below threshold eps
%   - only one cluster remains
%   - predefined number of clusters is reached
%
%for iter=1:maxIter
iter = 1;
runGravity = 1;
while runGravity
    
    % save iteration state
    if saveItersState == 1
        storedIters.data(iter,:) = {nClust, G, dataGrav, C, neighC, mass, valid};
    end
    
    % Stop criterion based on iterations
    % If maxIter is 0 ignore iteration check. Algorithms will stop due to
    % other criteria (reached predefined number of clusters, epsilon)
    if iter >= maxIter && maxIter ~= 0
        runGravity = 0;
    end

    % obviously: if number of clusters equals 1, then finish the job
    if nClust==1
        break;
    end
    
    % is the predefined number of clusters reached?
    if retVal.kPredefReached
        if showGrav
            fprintf(['Exiting - predefined number of clusters reached! k=',num2str(nClust),'\n']);
        end
        break;
    end
    
    % flag, which signalizes the convergence - the case, when movements of
    % particles are under threshold epsilon
    convergedGlobal=1;
    
    if(showGrav > 0)
        disp(['Iteration: ',num2str(iter),', G=',num2str(G),', clusters: ',num2str(nClust)]);
    end
    
    % random order of visiting particles
    visitOrder=valid(randperm(length(valid)));
    
    for i=visitOrder
        % if i-th particle doesn't exist, we skip it
        if ~isempty(C{i})
            
            % we select another particle k, neighbours of i are prefered in
            % accordance to value pOut (the smaller, the more neighbours are prefered)
            
            % list of neighbours of i-th particle
            ni=neighC{i};
            
            % if true, a neighbour of i is chosen (at random among all neighbours)
            if (~isempty(ni)) && (rand < 1-pOut)
                rpi=randperm(length(ni));
                % the chosen one
                k=ni(rpi(1));
                
                % with probability pOut, a random particle is chosen as k
            else
                choices=valid;
                choices(valid==i)=[];
                rpi=randperm(nClust-1);
                k=choices(rpi(1));
            end
            
            % Calculation of a new position of particles i and k.
            % Employ different strategies with respect to distance and
            % dist_mode.
            
            switch distance
                case 'sqEuclidean'
                    [dataGrav(i,:), dataGrav(k,:), converged, d_new] = gravityEuclidean(G,dataGrav(i,:),dataGrav(k,:),mass(i),mass(k),isMass,epsStop);
                    
                case 'correlation'
                    switch som_mode
                        case 'euc'
                            [dataGrav(i,:), dataGrav(k,:), converged, d_new] = gravityCorrelation(G,dataGrav(i,:),dataGrav(k,:),mass(i),mass(k),isMass,epsStop,'centered',grav_mode);
                            
                        case 'dotprod'
                            [dataGrav(i,:), dataGrav(k,:), converged, d_new] = gravityDotprodc(G,dataGrav(i,:),dataGrav(k,:),mass(i),mass(k),isMass,epsStop,grav_mode);
                    end
                case 'ucorrelation'
                    switch som_mode
                        case 'euc'
                            [dataGrav(i,:), dataGrav(k,:), converged, d_new] = gravityCorrelation(G,dataGrav(i,:),dataGrav(k,:),mass(i),mass(k),isMass,epsStop,'uncentered',grav_mode);
                        case 'dotprod'
                            [dataGrav(i,:), dataGrav(k,:), converged, d_new] = gravityDotprod(G,dataGrav(i,:),dataGrav(k,:),mass(i),mass(k),isMass,epsStop,grav_mode);
                    end
                    
                otherwise
                    error(['Wrong distance: ',distance]);
                    
            end
            if ~converged
                convergedGlobal = 0;
            end
            
            % check whether moved particles are close enough for merging
            % If d_new < 0 --> ignore, there was repulsion on work
            if(d_new < alfa)
                % update reference in C: add contents of cell k to cell i
                C{i} = [C{i} C{k}];
                
                % update mass; it is only relevant if flag isMass is set on 1 or 2
                mass(i)=mass(i)+mass(k);
                
                % clear the k particle cell
                C{k}=[];
                valid(valid==k)=[];
                
                % there is one particle less, thus one cluster
                % representative less
                nClust=nClust-1;
                
                
                % update neighbourhood info
                % k is not a neighbour of i anymore
                new=neighC{i};
                new(new==k)=[];
                neighC{i}=new;
                
                % all the neighbours of k are now neighbours of i
                nk=neighC{k};
                
                for nInd=1:length(nk)
                    new=neighC{nk(nInd)};
                    if ~ismember(i,new)
                        new(new==k)=i;
                    else
                        new(new==k)=[];
                    end
                    neighC{nk(nInd)}=new;
                    
                    if ~ismember(nk(nInd),neighC{i}) && nk(nInd) ~=i
                        neighC{i}=[neighC{i}, nk(nInd)];
                    end
                end
                
                %k has no neighbours
                neighC{k}=[];
                
            end
            
            % showtime - display each move of particle
            if(showGrav==1)
                
                pause(0.005);
                plot(dataStartOK(:,1),dataStartOK(:,2),'bo');
                axis('equal')
                hold on;
                
                if dimStart > 2
                    dataOK=(suns_pca_vec(:,1:2)' * dataGrav')';
                else
                    dataOK=dataGrav;
                end
                
                for j=1:nClust
                    
                    plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'k.','markersize',10);
                    plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'r.');
                    
                    for sosedI=1:length(neighC{C{valid(j)}(1)})
                        sourceIData = C{valid(j)}(1);
                        neighIData = neighC{C{valid(j)}(1)}(sosedI);
                        
                        if sourceIData < neighIData
                            plot(	[	dataOK(sourceIData,1); dataOK(neighIData,1)], ...
                                [	dataOK(sourceIData,2);dataOK(neighIData,2)],'k-');
                            
                        end
                    end
                    
                    title(['Iteration: ', num2str(iter),', G=',num2str(G)]);
                    text(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),num2str(mass(C{valid(j)}(1))),'verticalAlignment','bottom');
                end
                
                hold off;
                
                if dimStart <=2
                    switch distance
                        case 'sqEuclidean'
                            axis([0,1,0,1]);
                        case 'dotprod'
                            axis([-1,1,-1,1]);
                    end
                end
            end
            % is the predefined number of clusters reached?
            if nClust<=kPredef
                retVal.kPredefReached=1;
                break;
            end
        end
    end
    
    % if every move was is under threshold value, we break the loop
    if convergedGlobal
        if showGrav
            fprintf('Algorithm has converged!\n');
        end
        break;
        
    end
    % decrease G
    G = (1-dG)*G;
    
    % increment interation counter
    iter = iter + 1;
    
end
%MAIN LOOP - END

%------- RECONSTRUCTION OF LABELS -----------------------------------------
% final labels
target = zeros(size(dataStart,1),1);

% reconstruction of labels for data points
for c=1:nClust
    tmp=C{valid(c)};
    sL=winners.labels(tmp);
    for s=1:length(sL)
        target(sL(s)==sLabels)=c;
    end
end

% stop measuring time of gravitational algorithm
tG = toc(tidID);

% showtime - display final result
if(showGrav > 0)
    if(showGrav==2)
        
        plot(dataStartOK(:,1),dataStartOK(:,2),'ro');
        hold on;
        
        if dimStart > 2
            dataOK=(suns_pca_vec(:,1:2)' * dataGrav')';
        else
            dataOK=dataGrav;
        end
        
        for j=1:nClust
            plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'.');
        end
        
        hold off;
        
        if dimStart <=2
            axis([0,1,0,1]);
        end
        axis('square');
        title('Final positions of clusters representatives');
    end
       
    options.colorMode='mixed';
    options.title='Final result of gSOM';
    gSOM_scatterPlot(dataStart,target,nClust,options);
end



%------- RETURN VALUES ----------------------------------------------------
retVal.target = target;
retVal.nClust = nClust;
retVal.iter = iter;
retVal.time = [tSOM, tG, tSOM+tG];
retVal.SOMprop = winners.SOMprop;
retVal.SOMprop.bmus = winners.coords;
retVal.SOMtrain = winners.SOMtrain;
if saveItersState == 1
    retVal.SOMprop.sMap = sMap;
   retVal.storedIters = storedIters; 
end

%------- UNLOAD SOM TOOLBOX -----------------------------------------------
%remove SOM Toolbox from MATLAB path
rmpath('..\..\misc\somtoolbox_pplk');

end




function [winners,labels,sMap] = createSOM(data,msize,lattice,show,isMass,distance,advanced)
%[winners,labels]=createSOM(data,msize,shape,show,isMass,distance)
%
% gSOM - Level One - representing data with SOM
%
% INPUTS
%   data:		nD x d matrix of data, where nD is the number of data points
%				and d the dimensionality of data points.
%   msize:		size of map M=a x b (1x2 vector of a and b), e.g. [17,12];
%				if empty [], heuristic is used M=ceil(5*sqrt(nD)), where nD is the
%				number of data points;
%				map size can also be determined with provided scale factor S.
%				Toolbox recognizes this by negative sign, e.g. use msize=[-0.75] to
%				apply size of 0.75*5*sqrt(nD). Default is -1.
%	lattice:	shape of the SOM grid ('rect' for rectangular - default, 'hexa' for
%				hexagonal)
%	show:		if 1, trained SOM is displayed in a figure
%	isMass:		[experimental!] if 1, mass of neurons is computed as number of data points that
%				belong to certain winning neuron. Default is 0.
%	distance:	distance metric: 'sqEuclidean', 'correlation', 'ucorrelation'
%	advanced:	additional settings for training, given as name-value cell array.
%               options:
%					- 'algorithm'	('seq' | 'batch')
%					- 'radius_ini'	(numeric)
%					- 'radius_fin'	(numeric)
%					- 'alpha_ini'	(numeric)
%					- 'alpha_type'	('linear' | 'inv' | 'power')
%					- 'initalg'		('randinit' | 'lininit')
%					- 'trainlen'	('short' | 'long' | [numRoughPhase, numFinePhase])
%					- 'sampleOrder'	('ordered' | 'random')
%					- 'neigh'		('bubble' | 'gaussian' | 'cutgauss' | 'ep')
%					- 'normalize'	(0 | 1)
%
% OUTPUTS
%   winners :   stucture with fields
%					.SOMprop	: properties of SOM
%					.SOMtrain	: train info about SOM
%                   .labels		: labels of SOM units (to link with data points)
%                   .coords		: coordinates of SOM units
%                   .neigh		: pairs of neighbours (a,b), where a < b
%                   .mass		: mass (number of data points) of each winning neuron
%					.pca_vec	: if data are of dimensionality > 2, there
%								would be eigen vectors computed for winners
%					.dataPca	: PCA of data, if dim. > 2
%
%   labels :   labels of winning neurons - each data point has its
%              coresponding winning neuron.
%
%--------------------------------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
%==========================================================================

if ~exist('isMass','var')
    isMass=0;
end

if ~exist('show','var')
    show=0;
end

if ~exist('advanced','var')
    advanced=[];
else
    i=1;
    while i<=length(advanced)
        if ischar(advanced{i})
            switch advanced{i}
                case 'normalize'
                    i=i+1; normalize = advanced{i};
                otherwise
                    i=i+1;
            end
        end
        i=i+1;
    end
end
if ~exist('normalize','var')
    normalize = 0;
end

[nD,d]=size(data);

D=som_data_struct(data);
D.isNormalized = normalize;

% SOM Toolbox function som_make default is batch mode of the SOM training
if any(strcmpi(lattice,{'rect','hexa'}))
    if isempty(msize)
        sMap = som_make(D,'lattice', lattice,'distfun',distance,'tracking',show,'advanced',advanced);
    else
        sMap = som_make(D,'msize',msize,'lattice', lattice,'distfun',distance,'tracking',show,'advanced',advanced);
    end
end

winners.SOMprop=sMap.topol;
winners.SOMprop.codebook=sMap.codebook;
winners.SOMtrain=sMap.trainhist;
winners.advancedParams=advanced;

sTrain=sMap.trainhist(length(sMap.trainhist));

% Find winning neurons (BMUs) for data points
if isfield(sTrain,'bmus') && ~isempty(sTrain.bmus) && sum(sTrain.bmus) > 0
    labels = sTrain.bmus;
else
    labels=som_bmus(sMap, D);
end
%labels=som_bmus(sMap, D);

%determine unique BMUs -> these are winning neurons labels
labelsUniq=unique(labels);

%how many data points are represented by certain winning neuron?
mass=hist(labels,labelsUniq);


%save labels and positions of winning neurons
winners.labels=labelsUniq;
winners.coords = sMap.codebook(labelsUniq,:);

if(isMass)
    winners.mass=mass;
else
    winners.mass=ones(1,length(mass));
end


%compute neighbours
N = som_connection(sMap);

%create matrix Npairs x 2 for storing neighbourhood info
pairs = nchoose2(labelsUniq);

neigh=[];

%check whether certain pair is neighbouring
for ind=1:size(pairs,1)
    if(N(pairs(ind,1),pairs(ind,2))==1)
        neigh=[neigh; pairs(ind,:)];
    end
    
end

% create cell struct - each sun has a list of its neighbours (indices in
% labelsUniq)
neighC = cell(1,length(labelsUniq));
for s=1:length(labelsUniq)
    if isempty(neigh)
        x=[];
        y=[];
    else
        x=find(neigh(:,1)==labelsUniq(s));
        y=find(neigh(:,2)==labelsUniq(s));
    end
    
    if ~isempty(x)
        nX=neigh(x,2);
        for nXi=1:length(nX)
            neighC{s}=[neighC{s}, find(labelsUniq==nX(nXi))];
        end
    end
    if ~isempty(y)
        nY=neigh(y,1);
        for nYi=1:length(nY)
            neighC{s}=[neighC{s}, find(labelsUniq==nY(nYi))];
        end
    end
end

winners.neigh=neighC;

%------- PLOT -------------------------------------------------------------
% eigen vectors of winning neurons, required for plotting
winners.pca_vec=[];
% PCA mapping of the data
winners.dataPca=[];

if(show)
    % number of winner codebooks
    numCB=size(sMap.codebook,1);
    
    if d > 2
        %PCA mapping of neurons
        neur=[sMap.codebook ; winners.coords];
        %first two principal components are already computed
        new_CB=(D.pca_vec(:,1:2)' * neur')';
        
        CB=new_CB(1:numCB,:);
        suns_new=new_CB(numCB+1:end,:);
        
        %PCA mapping of data
        new_data=(D.pca_vec(:,1:2)' * D.data')';
        
        winners.pca_vec=D.pca_vec;
        winners.dataPca=new_data;
        
    else
        CB=sMap.codebook;
        suns_new=winners.coords;
        new_data=D.data;
        
    end
    
    fig=figure();
    options.fig=fig;
    gSOM_scatterPlot(new_data,labels,numCB,options);
    hold on;
    
    title('Level One - SOM');
    som_grid(sMap, 'Coord', [CB(:,1),CB(:,2)]);
    axis('equal')
    
    
    %plot(new_data(:,1),new_data(:,2),'bo');
    
    for i=1:length(winners.mass)
        plot(suns_new(i,1),suns_new(i,2),'r.', 'markersize', winners.mass(i)/max(winners.mass)*20)
        text(suns_new(i,1),suns_new(i,2)+0.01,num2str(winners.mass(i)),'verticalAlignment','bottom');
    end
    hold off;
    
end
end

function [x_new, y_new, converged, d_new]=gravityEuclidean(G,x,y,mass_x,mass_y,isMass,epsStop)
% R = y-x   (vector pointing towards y)
R = y - x;
dist_x2y = norm(R);

% E is unit vector of R
E = R ./ dist_x2y;

% force is proportional to mass
if(isMass==1)
    F=(E.*G*mass_x*mass_y)./(dist_x2y^2);
    
% inverse - greater mass, lower gravity force
elseif(isMass==2)
    F=(E.*G)./(dist_x2y^2 *mass_x*mass_y);
    
% default is to neglect the mass; proves to be OK
else
    F=(E.*G)./(dist_x2y^2);
end

% the move of each particle (i and k) is limited to half of distance
% between the pair
move = min( dist_x2y /2 , norm(F) );

% movement of each particle has to be less than epsilon to stop the algorithm
if move > epsStop
    converged=0;
else
    converged=1;
end

% perform the move of selected particles
x_new = x + move.*E;
y_new = y - move.*E;

% new distance between them
d_new = norm(x_new-y_new);

end

function [x_new, y_new, converged, d_new]=gravityCorrelation(G,x,y,mass_x,mass_y,isMass,epsStop,distance,grav_mode)
% distance = {'centered', 'uncentered'
% grav_mode = {'AR' | 'semi' | 'pure'}

% R = y-x   (vector pointing towards y)
R = y - x;

[r, dist_x2y] = pearson_corr(x,y,distance,0); % vary last two arguments for other correlation functions

if strcmp(grav_mode,'AR')
    dist_x2y=norm(R); %euclidean distance
end

% E is an unit vector of R
E = R ./ dist_x2y;

% force is proportional to mass
if(isMass==1)
    F=(E.*G*mass_x*mass_y)./(dist_x2y^2);
    
    % inverse - greater mass, lower gravity force
elseif(isMass==2)
    F=(E.*G)./(dist_x2y^2 *mass_x*mass_y);
    
    %default is to neglect the mass; proves to be OK
else
    F=(E.*G)./(dist_x2y^2);
end

% the move of each particle (i and k) is limited to half of distance between
% the pair

switch grav_mode
    case 'AR'
        move = r * min( dist_x2y /2 , norm(F) ); % r is from [-1,1]
        
        % perform the move of selected particles
        x_new = x + move.*E;
        y_new = y - move.*E;
        
        % new distance between them
        d_new = norm(x_new-y_new);
        
    case 'semi'
        move = min( dist_x2y /2 , norm(F) );
        
        % perform the move of selected particles
        x_new = x + move.*E;
        y_new = y - move.*E;
        
        % new distance between them
        [r_dummy, d_new] = pearson_corr(x_new,y_new,distance,0);
        
    case {'pure' , 'default'}
        [r_move,d_move] = pearson_corr(x,x+F,distance,0);
        move = min( dist_x2y /2, d_move);
        
        % perform the move of selected particles
        x_new = x + move.*E;
        y_new = y - move.*E;
        
        % new distance between them
        [r_dummy, d_new] = pearson_corr(x_new,y_new,distance,0);
        
    otherwise
        error(['Wrong mode: ',mode]);
end

% movement of each particle has to be less than epsilon to stop the algorithm
if move > epsStop
    converged=0;
else
    converged=1;
end

end

function [x_new, y_new, converged, d_new] = gravityDotprod(G,x,y,mass_x,mass_y,isMass,epsStop,mode)

% mode: {'linear' | 'sqLinear' | 'symmetric' | 'sqSymmetric'}

% We assume that points x and y are already normalized (in SOM training)
% and the vectors that represents them have unit length. So, they are
% situated somewhere on a unit sphere and distance between x and y is
% actually the angle fi between them - dot product (or uncentered correlation)
% of normalized vectors is nothing more than cos(fi).
%
% According to the angle fi and G (optionally mass), we rotate the vector x
% to/from y and y to/from x to simulate force of gravity or repulsion,
% respectively.
%		Gravity:	0  < cos(fi) < 1	-90 < fi < 90
%		Repulsion:	-1 < cos(fi) < 0	 90 < fi < 270

% r = cos(fi) = (x * y') / (||x|| * ||y||)
% assuming: ||x|| = ||y|| = 1  => cos(fi) = x * y'

% difference of vectors
R = y - x;

r = x*y'; % r = cos(fi)

% direction of movement (gravity/repulsion)
%	r > 0 -> dir = +1 -> gravity
%	r < 0 -> dir = -1 -> repulsion
%	r = 0 -> dir =  0 -> stay still
% RELEVANT ONLY FOR mode='symmetric' or 'sqSymmetric'

dir=1;


% distance between x and y -> norm of R
% it can be:
%		2 * (1 -  r  )  --> points far away (negative r) will move
%							really slowly away from each other
%		or
%		2 * (1 - |r| )  --> more than r is negative, the stronger repulsion
%							between x and y gets.
%
% try:
% fi=0:pi/20:pi; r=cos(fi); plot(r,2*(1-r),'-.'); hold on; plot(r,1*(1-abs(r)));
switch mode
    case {'linA', 'default'}
        dist_x2y = 2*(1-r);
        
    case 'sqA'
        dist_x2y = (2*(1-r))^2;
        
    case 'linAR'
        dir=sign(r);
        dist_x2y = 2*(1-dir*r);
        
    case 'sqAR'
        dir=sign(r);
        dist_x2y = (2*(1-dir*r))^2;
        
    otherwise
        error(['Wrong mode: ',mode]);
end

% force is proportional to mass
if(isMass==1)
    F = (G * mass_x * mass_y) / (dist_x2y + eps);
    
    % inverse - greater mass, lower gravity force
elseif(isMass==2)
    F = G /(dist_x2y * mass_x * mass_y);
    
    %default is to neglect the mass; proves to be OK
else
    F = G / (dist_x2y + eps);
end

move = min(0.5, F);

% movement of each particle has to be less than epsilon to stop the algorithm
if move > epsStop
    converged=0;
else
    converged=1;
end

x_new = x + dir*move*R;
y_new = y - dir*move*R;

%normalize x and y
x_new = x_new ./ norm(x_new);
y_new = y_new ./ norm(y_new);

% new distance between x and y
r = x_new * y_new';

if dir == -1
    d_new = Inf; % to prevent merging of two anti-correlated points
    
else
    switch mode
        case {'linA', 'default'}
            d_new = 2*(1-r);
        case 'sqA'
            d_new = (2*(1-r))^2;
            
        case 'linAR'
            d_new = 2*(1-abs(r)); % if d_new < 0 -> repulsion, DO NOT MERGE
        case 'sqAR'
            d_new = (2*(1-abs(r)))^2;  % if d_new < 0 -> repulsion, DO NOT MERGE
    end
end

end

function [x_new, y_new, converged, d_new] = gravityDotprodc(G,x,y,mass_x,mass_y,isMass,epsStop,mode)

% mode: {'linA' | 'sqA' | 'linAR' | 'sqAR'}

% We assume that points x and y are already normalized (in SOM training)
% and the vectors that represents them have unit length. So, they are
% situated somewhere on a unit sphere and distance between x and y is
% actually the angle fi between them - dot product (or uncentered correlation)
% of normalized vectors is nothing more than cos(fi).
%
% According to the angle fi and G (optionally mass), we rotate the vector x
% to/from y and y to/from x to simulate force of gravity or repulsion,
% respectively.
%		Gravity:	0  < cos(fi) < 1	-90 < fi < 90
%		Repulsion:	-1 < cos(fi) < 0	 90 < fi < 270

% r = cos(fi) = ( (x-mean(x)) * (y-mean(y))') / (||x-mean(x)|| * ||y-mean(y)||)
% assuming: ||x|| = ||y|| = 1  => cos(fi) = x * y'

% difference of vectors
R = y - x;

r = x*y'; % r = cos(fi)

% direction of movement (gravity/repulsion)
%	r > 0 -> dir = +1 -> gravity
%	r < 0 -> dir = -1 -> repulsion
%	r = 0 -> dir =  0 -> stay still
% RELEVANT ONLY FOR mode='symmetric' or 'sqSymmetric'

dir=1;


% distance between x and y -> norm of R
% it can be:
%		2 * (1 -  r  )  --> points far away (negative r) will move
%							really slowly away from each other
%		or
%		2 * (1 - |r| )  --> more than r is negative, the stronger repulsion
%							between x and y gets.
%
% try:
% fi=0:pi/20:pi; r=cos(fi); plot(r,2*(1-r),'-.'); hold on; plot(r,1*(1-abs(r)));
switch mode
    case {'linA', 'default'}
        dist_x2y = 2*(1-r);
        
    case 'sqA'
        dist_x2y = (2*(1-r))^2;
        
    case 'linAR'
        dir=sign(r);
        dist_x2y = 2*(1-dir*r);
        
    case 'sqAR'
        dir=sign(r);
        dist_x2y = (2*(1-dir*r))^2;
        
    otherwise
        error(['Wrong mode: ',mode]);
end

% force is proportional to mass
if(isMass==1)
    F = (G * mass_x * mass_y) / (dist_x2y + eps);
    
    % inverse - greater mass, lower gravity force
elseif(isMass==2)
    F = G /(dist_x2y * mass_x * mass_y);
    
    %default is to neglect the mass; proves to be OK
else
    F = G / (dist_x2y + eps);
end

move = min(0.5, F);

% movement of each particle has to be less than epsilon to stop the algorithm
if move > epsStop
    converged=0;
else
    converged=1;
end

x_new = x + dir*move*R;
y_new = y - dir*move*R;

%normalize x and y
x_new = (x_new-mean(x_new)) ./ norm((x_new-mean(x_new)));
y_new = (y_new-mean(y_new)) ./ norm((y_new-mean(y_new)));

% new distance between x and y
r = x_new * y_new';

if dir == -1
    d_new = Inf; % to prevent merging of two anti-correlated points
    
else
    switch mode
        case {'linA', 'default'}
            d_new = 2*(1-r);
        case 'sqA'
            d_new = (2*(1-r))^2;
            
        case 'linAR'
            d_new = 2*(1-abs(r)); % if d_new < 0 -> repulsion, DO NOT MERGE
        case 'sqAR'
            d_new = (2*(1-abs(r)))^2;  % if d_new < 0 -> repulsion, DO NOT MERGE
    end
end

end