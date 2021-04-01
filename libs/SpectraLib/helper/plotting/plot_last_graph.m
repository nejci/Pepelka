function plot_last_graph(plot_stdev,outfile) 
% function plot_last_graph(plot_stdev,outfile) 
% Plots the last graph stored in outfile. ('default,
% tmp/last_graph_plotted

  if nargin<1
      plot_stdev=0;
  end
  if nargin<2
	outfile='tmp/last_graph_plotted';
  end
  
  load(outfile,'x','y','y_stdev','CS','DATA_FILE','metric')
  if ~plot_stdev
      y_stdev=[];
  end
  
  plot_graph_generic(x,y,y_stdev,CS,strcat(DATA_FILE,' METRIC=',metric));
  
