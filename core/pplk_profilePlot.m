function [fig]=pplk_profilePlot(data,labels,K,varargin)
% [fig]=pplk_profilePlot(data,labels,K,varargin)
%
% INPUTS		
%   data
%       A N-by-D matrix of input data.
%
%   labels
%       A vector of data labels - clustering.
%
%   options			
%       Options structure or name-value list.
%       It's fields are:
%           - figTitle
%           - axisLabels a 1-by-2 cell of strings
%           - axisTicks a 1-by-2 cell of vectors
%           - axisTicksLabels a 1-by-2 cell of string 1-by-D cells
%           - view ['overview'|'scatterPlot'|'details'|'PCA']
%           - colorMode ['color'|'pattern'|'mixed'|'mixedClear'] 
%             (only for scatterPlot)
%           - normalize - (scatterPlot) normalize on interval [0,1]
%           - fig - force this figure handle
%           - annotations a n-by-1 cell of strings, where n is the number 
%             of data points
%           - centroids [0 | 1] - plot cluster centroid or not 
%
%
% OUTPUTS	
%   fig
%       - A handle on figure.
%       - An array of handles on figures.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

callDir=chdir(pplk_homeDir());

options=varargin;

%DEFAULTS
figTitle='';
axisLabels={'',''};
axisTicks={[],[]};
axisXTicksLabels={''};
colorMode = 'color';
view='overview';
normalize=0;
annotations={};
centroids = 1;

% constants
COLOR_GENES=[0.8, 0.8, 0.8];
COLOR_CENTROIDS = 'magenta';
WIDTH_CENTROIDS = 2;
COLOR_ERRORBAR = 'black';
WIDTH_ERRORBAR = 2;

if ~exist('data','var')
	error('Too few arguments - data missing!');
end

[N,D]=size(data);
Yval_min = min(min(data));
Yval_min = sign(Yval_min)*abs(ceil((-1)*Yval_min));

Yval_max = max(max(data));
Yval_max = sign(Yval_max)*abs(ceil(Yval_max));

if ~exist('labels','var') || isempty(labels)
	labels = ones(N,1); % one cluster of all data points
	K=1;
end

if ~exist('K','var') || isempty(K)
	%number of clusters
	K = 1;
end

if(exist('options','var'))
	if isstruct(options)
		if isfield(options,'title')
			figTitle=options.title;
		end
		if isfield(options,'axisLabels')
			axisLabels=options.axisLabels;
		end
		if isfield(options,'axisTicks')
			axisTicks=options.axisTicks;
		end
		if isfield(options,'axisXTicksLabels')
			axisXTicksLabels=options.axisXTicksLabels;
		end
		if isfield(options,'view')
			view=options.view;
		end
		if isfield(options,'colorMode')
			colorMode=options.colorMode;
		end
		if isfield(options,'normalize')
			normalize=options.normalize;
		end
		if isfield(options,'fig')
			fig = options.fig;
		end
		if isfield(options,'annotations')
			annotations=options.annotations;
		end
		if isfield(options,'centroids')
			centroids=options.centroids;
		end
	else
		
		if ~isempty(options)
			i=1;
			while i<=length(options)
				if ischar(options{i})
					switch options{i}
						case 'title', i=i+1; figTitle = options{i};
						case 'axisLabels', i=i+1; axisLabels = options{i};
						case 'axisTicks', i=i+1; axisTicks = options{i};
						case 'axisXTicksLabels', i=i+1; axisXTicksLabels = options{i};
						case 'view', i=i+1; view = options{i};
						case 'colorMode', i=i+1; colorMode = options{i};
						case 'normalize', i=i+1; normalize = options{i};
						case 'fig', i=i+1; fig = options{i};
						case 'annotations', i=i+1; annotations = options{i};
						case 'centroids', i=i+1; centroids = options{i};
						otherwise
							warning(['Option value ''',options{i},''' is not supported!']);
							i=i+1;
					end
				end
				i=i+1;
			end
		end
		
	end
		
end

if ~exist('fig','var') || isempty(fig)
	fig=figure();
else
	figure(fig);
end
%-------------------------------------------------------------------------


switch view
	
	case 'scatterPlot'
		% pass the arguments to the pplk_scatterPlot
		optTmp.title = figTitle;
		optTmp.axisLabels = axisLabels;
		optTmp.axisTicks = axisTicks;
		optTmp.colorMode = colorMode;
		optTmp.normalize = normalize;
		optTmp.annotations = annotations;
		if exist('fig','var'), optTmp.fig = fig; end
		[fig]=pplk_scatterPlot(data,labels,K,optTmp);
		
	case 'overview'
		gridWidth = ceil(sqrt(K));
		gridHeight = ceil(K/gridWidth);
		
		if gridHeight > gridWidth
			tmp=gridHeight;
			gridHeight = gridWidth;
			gridWidth = tmp;
		end
		
		if(K<5)
			colors=[0 0 153/255; 0 153/255 1; 153/255 1 102/255; 1 102/255 0];
			
		else
			colors=colormap(jet(K*10));
			colors=colors(1:10:end,:);
		end
		
		fig2=figure(); plot(0,0); h2=get(fig2,'CurrentAxes');
		fig=[fig, fig2];
		
		for k=1:K
			figure(fig(1));
			h=subplot(gridHeight,gridWidth,k);
			%set(h,'Units','pixels');
			%p=get(h,'pos');
			%disp(p);
			%set(h,'pos',p+[-10, 0, 20, 20]);
			%set(h,'Units','Normalized');
			
			cluster_i = data(labels==k,:);
			plot(cluster_i','LineStyle','-','Color',COLOR_GENES);
			text(0.98, 0.96,[num2str(size(cluster_i,1)),' genes'],'HorizontalAlignment','right','FontName','FixedWidth','Units','Normalized');
			
			
			% calculate mean and std for centroids
			if centroids
				hold on;
				centr_mean = mean(cluster_i,1);
				centr_std = std(cluster_i,[],1);
				errorbar(centr_mean,centr_std,'Color',COLOR_ERRORBAR,'LineWidth',WIDTH_ERRORBAR);
				plot(centr_mean,'LineStyle','-','Marker','.','Color',colors(k,:),'LineWidth',WIDTH_CENTROIDS);
				text(0.98, 0.04,['STD=',num2str(sum(centr_std))],'HorizontalAlignment','right','FontName','FixedWidth','Units','Normalized');
				hold off;
			end
			
			title(h,['cluster ',num2str(k)]);
			
			if isempty(figTitle)
				titleHnd=get(h2,'title');
				set(titleHnd,'String','Overview');
			else 
				titleHnd=get(h2,'title');
				set(titleHnd,'String',[figTitle,' (overview)']);
			end
			
			% X axis
			if ~isempty(axisTicks{1})
				set(h, 'XTick',axisTicks{1});
				set(h2,'XTick',axisTicks{1});
			else
				set(h, 'XTick',1:1:D);
				set(h, 'XLim',[1,D]);
				
				set(h2, 'XTick',1:1:D);
				set(h2, 'XLim',[1,D]);
			end
			
			if ~isempty(axisXTicksLabels{1})
				set(h, 'XTickLabel',axisXTicksLabels);
				set(h2, 'XTickLabel',axisXTicksLabels);
			else
				set(h, 'XTickLabel','');
				set(h2, 'XTickLabel','');
			end
			
			% Y axis
			if ~isempty(axisTicks{2})
				set(h, 'YTick',axisTicks{2});
				set(h2, 'YTick',axisTicks{2});
			else
				set(h, 'YLim',[Yval_min,Yval_max]);
				set(h, 'YMinorTick','on');
				
				set(h2, 'YLim',[Yval_min,Yval_max]);
				set(h2, 'YMinorTick','on');
			end
			
			figure(fig(2));
			hold on;
			plot(cluster_i','LineStyle','-','Color',colors(k,:));
		end
		
		figure(fig(2));
		hold off;
					
		
	case 'details'
		
	case {'PCA','pca'}
		
		%[pc,newData]=princomp(data,'econ');		
		%plot(newData(:,1),newData(:,2),'b.');
		%biplot(pc(1:2,1:2),'scores',newData(:,1:2),'ObsLabels',annotations);
		% 		if isempty(figTitle)
		% 			title('PCA');
		% 		else
		% 			title([figTitle,' (PCA)']);
		% 		end
		
		close(fig);
		if isempty(annotations)
			mapcaplot(data);
		else
			mapcaplot(data, annotations);
		end


	otherwise
		error(['Wrong view mode: ''',view,'''']);
	
end

chdir(callDir);