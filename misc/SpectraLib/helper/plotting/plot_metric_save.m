function plot_metric_save(DATA_FILE,CS,KRANGE,iterations,metric,plot_stdev,outdir)
% function plot_metric_save(DATA_FILE,CS,KRANGE,iterations,metric,plot_stdev,outdir)
% 
% Plots the graph for given data file for CS algo and KRANGE and
% NUMBER_OF_RUNS and metric 
SIGMA=-1; 

clear x y y_stdev

if nargin < 6
		plot_stdev=1;
end

if nargin < 7
  outdir='/scratch/deepak/spectral/newout/';
else
  outdir=strcat(outdir,'/');
end


mink=min(KRANGE); 
for algoiter=1:length(CS)
  for k=KRANGE
	stats=stats_single_experiment(DATA_FILE,CS{algoiter},k,SIGMA,iterations,outdir);
	if metric=='ce'
	  curr_stat=stats(1);
	  curr_stdev=stats(2);
      metric_name='Clustering Error'; 
	end
	if metric=='vi'
	  curr_stat=stats(3);
	  curr_stdev=stats(4);
      metric_name='Variation of Information';
	end
	if metric=='ai'
	  curr_stat=stats(5);
	  curr_stdev=stats(6);
      metric_name='Adjusted Rand Index';
	end	
	if metric=='wi'
	  curr_stat=stats(7);
	  curr_stdev=stats(8);
      metric_name='Wallas Index (One sided)'; 
	end	
	y(algoiter,k-mink+1)=curr_stat; 
	y_stdev(algoiter,k-mink+1)=curr_stdev; 
  end
end
x=KRANGE; 

outfile='tmp/last_graph_plotted'; 
  
save(outfile,'x','y','y_stdev','CS','DATA_FILE','metric')
if plot_stdev
		y_stdev=[];
end

plot_graph_generic(x,y,y_stdev,CS,strcat(DATA_FILE,' '));
% plot_graph_generic(x,y,y_stdev,CS,strcat(DATA_FILE,' [ ',metric_name,' ]'));

xlabel('Number of Clusters (K)');
ylabel(metric_name); 
