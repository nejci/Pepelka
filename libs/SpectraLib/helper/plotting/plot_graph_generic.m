function plot_graph_generic(x,y, y_stdev ,legend_string, title_string,logx,logy,color_vector,noclf)
%function plot_graph_generic(x,y, y_stdev ,legend_string,title_string,logx,logy,color_vector,noclf)
% 
% Generic function to plot a graph. Each row of y represents a line to
% plot 
  
  if nargin < 3 
	y_stdev =''; 
  end
  if nargin <4 
	legend_string=''; 
  end
  if nargin < 5
	title_string=''; 
  end
  if nargin < 7
	logx=0; 
	logy=0; 
  end
  if nargin < 8 
	color_vector=''; 
  end
  if nargin < 9
	noclf=0; 
  end
  
  if length(color_vector)==0
	global DEFAULT_COLOR_VECTOR; 
	color_vector=DEFAULT_COLOR_VECTOR; 
  end
  
  number_of_points=length(x); 
  number_of_lines=size(y,1) ;% number of rows. 
  
  if (~noclf)
	clf; 
  end
  
  
  for i=1:number_of_lines 
	if i==2 
	  hold on 
	end
	
	if logx
	  if logy 
		loglog(x,y(i,:),color_vector(i,:) ); 
	  else	  
		semilogx(x,y(i,:),color_vector(i,:) ); 
	  end
	else
	  if logy 
		semilogy(x,y(i,:),color_vector(i,:) ); 
	  else
		plot(x,y(i,:),color_vector(i,:),'LineWidth',2 );  % the default case 
        set(gca,'xtick',x) 
	  end	
	end	  
  end
  if length(y_stdev)>0
	for i=1:number_of_lines 
	  errorbar(x,y(i,:),y_stdev(i,:),color_vector(i,:) );  
	
	end
  end

  set(gca,'FontSize',15); 
  
  if length(legend_string)  ~= 0
	legend(strrep(legend_string,'_','-'),0); 
  end
  
  set(gca,'FontSize',17); 

  
  if length(title_string) ~=0
	title(strrep(title_string,'_','-')); 
  end
  
