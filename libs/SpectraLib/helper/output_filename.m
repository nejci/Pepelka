function filename=output_filename(dir,dataFilePrefix,cluster_algo,k,sigma,iteration_number)
% function filename=output_filename(dir,dataFilePrefix,cluster_algo,k,sigma,iteration_number)
% The output filename based on the various parameters. 
% NOTE: The directory must have the slash in it.   

  filename=strcat(dir,dataFilePrefix,'_',cluster_algo,'_K',sprintf('%d',k),'_',sigma_str(sigma),'_',sprintf('%d',iteration_number));
  
	
