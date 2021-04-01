function stats=stats_single_experiment(dataFile,cluster_algo,k,sigma,iterations,outdir) 
% function stat_single_experiment(dataFile,cluster_algo,k,iterations,sigma,outdir,plot_points)
%
% Function which reads the output of run_single_experiment and then
% returns the stats. 
  if nargin == 5
	outdir='newout/';
  end
  niter=length(iterations); 
  clear vi ce ari; 
  vi=zeros(niter,1); 
  ce=zeros(niter,1); 
  ari=zeros(niter,1); 
  wi=zeros(niter,1); 
  [S,true_ca]=read_from_data_file(dataFile); 
  for iteration_number=iterations 
	% this is what makes one iteration different from another. 
	
	outFileName=output_filename(outdir, dataFile, cluster_algo, k, sigma, iteration_number) ;
	load(outFileName, 'variation_of_information','clustering_error_true','adjRI','cluster_assignment');
	vi(iteration_number)=variation_of_information; 
	ce(iteration_number)=clustering_error_true; 
	ari(iteration_number)=adjRI; 
	wi(iteration_number)=wallas_index(true_ca,cluster_assignment); 
	
  end
  
  stats(1)=mean(ce); 
  stats(2)=std(ce); 
  stats(3)=mean(vi); 
  stats(4)=std(vi); 
  stats(5)=mean(ari); 
  stats(6)=std(ari); 
  stats(7)=mean(wi); 
  stats(8)=std(wi); 
  
  