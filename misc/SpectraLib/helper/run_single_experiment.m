function run_single_experiment(dataFile,cluster_algo_list,k_range,sigma,iterations,outdir,plot_points) 
% function run_single_experiment(dataFile,cluster_algo_list,k_range,sigma,iterations,outdir,plot_points) 
%
% Function which runs a single experiment and saves the result in an
% appropriate file. The S is assumed to be given in the file. 
% iterations can be a single number or a list in which case so
% many iterations are run  
% if the dataFile is a cell array then it treats the first name as the
% output file and the second as the similarity matrix Matrix it treats it
% like a similarity matrix and the first is the filename. 
  
  tic; 
  if nargin < 6
	outdir='/scratch/deepak/spectral/newout/';
  else
	outdir=strcat(outdir,'/'); 
  end
  if nargin < 7 
	plot_points=0;
  end
  SINGLE_VERBOSE=0; 
  if iscell(dataFile) 
	  dataFileName=dataFile{1}; 
	  similarity=dataFile{2}; 
  else
      dataFileName=dataFile;
	% read in the input data
	  [similarity,orig_cluster,points]=read_from_data_file(dataFile);
	  
	  
	% len_points=length(points)
	% calculate the similarity in case sigma > 0 
	  if (length(points) && sigma > 0)
	similarity=S_from_points(points,sigma,smoothing_constant,0); 
	  end
  end 
  niter=length(iterations); 
  
  for algoiter=1:length(cluster_algo_list)
	cluster_algo=cluster_algo_list{algoiter};
	for k=k_range
	  fprintf('Running %s (%d iterations) at %.2fs (algoiter=%d)\n',output_filename('',dataFileName,cluster_algo,k,sigma,nan),niter,toc,algoiter);
	  clear vi ce ari; 
	  vi=zeros(niter,1); 
	  ce=zeros(niter,1); 
	  ari=zeros(niter,1); 
	  for iteration_number=iterations 
		% this is what makes one iteration different from another. 
		rand('state',iteration_number); 
		
		output_str=''; 
		outFileName=output_filename(outdir, dataFile, cluster_algo, k,sigma, iteration_number); 
		if SINGLE_VERBOSE
		  fprintf('Writing to file=%s\n',outFileName);
		end   
		fp=fopen(strcat(outFileName,'.out'),'w');
		
		if fp==-1
		  if SINGLE_VERBOSE
			fprintf('Cannot open file=%s',filename);
		  end
		  dbstop; 
		end

		% need to write the following
		%
		% dataFile
		% cluster_algo
		% options to the cluster_algo 
		% classification by the cluster algo.
		% variation of information
		% clustering error
		% permutation for which clustering error is minimized.
		
		
		% SMOOTHING CONSTANT. = -1 in case no smoothing. 
		smoothing_constant=-1; 
		
		output_str=strcat(output_str,'\n',sprintf('dataFile=%s\n',dataFile));
		output_str=strcat(output_str,'\n',sprintf('cluster_algo=%s\n',cluster_algo));
		output_str=strcat(output_str,'\n',sprintf('k=%d\n',k));
		output_str=strcat(output_str,'\n',sprintf('sigma=%.3f\n',sigma));
		output_str=strcat(output_str,'\n',sprintf('smoothing_constant=%g\n',smoothing_constant));

		

		%   % run the clustering algo
		cluster_assignment=feval(cluster_algo,similarity,k);
		
		% write out the information to the file 
		output_str=strcat(output_str,'\n',sprintf('cluster_algo_options_begin\n'));  
		%	output_str=strcat(output_str,'\n',sprintf( '%s\n',param_string));
		output_str=strcat(output_str,'\n',sprintf('cluster_algo_options_end\n'));  
		output_str=strcat(output_str,'\n',sprintf('\n'));
		
		
		variation_of_information=-1;
		clustering_error_true=-1;
		% calculate/print the confusion matrix and the "distance" in case there
		% are any true clusters. 
		if (~isempty(orig_cluster))
		  %variation of information 
		  variation_of_information=compare_clusterings(orig_cluster, cluster_assignment);
		  vi(iteration_number)=variation_of_information; 
		  output_str=strcat(output_str,'\n',sprintf('variation_of_information=%-.4f\n',variation_of_information));
		  if SINGLE_VERBOSE
			fprintf('variation_of_information=%-.4f\n',variation_of_information);
		  end
		  %clustering error 
		  [clustering_error_true,permutation]=clustering_error(orig_cluster, cluster_assignment);
		  output_str=strcat(output_str,'\n',sprintf('clustering_error_true=%-.4f\n',clustering_error_true));
		  ce(iteration_number)=clustering_error_true; 
		  if SINGLE_VERBOSE
			fprintf('clustering_error_true=%-.4f\n',clustering_error_true);
		  end
		  output_str=strcat(output_str,'\n',sprintf('permutation=%s\n',sprintf(' %d',permutation)));
		  
		  % Adjusted Rand Index 
		  adjRI=adjusted_rand_index(orig_cluster, cluster_assignment); 
		  ari(iteration_number)=adjRI; 
		  output_str=strcat(output_str,'\n',sprintf('adjRI=%-.4f\n',adjRI));
		  if SINGLE_VERBOSE
			fprintf('Adjusted Rand Index=%-.4f\n',adjRI);
		  end
		else % print out the cluster indices. 
		  output_str=strcat(output_str,'\n',sprintf('\n\nERROR NO ORIGINAL CLUSTERING FOUND\n\n'));
		  %  for i=1:k
		  %	disp(fprintf('Cluster %d',i));
		  %	disp(find(cluster_assignment==i))
		  %  end
		end

		% write down the classification in any case.
		
		output_str=strcat(output_str,'\n',sprintf('cluster_assignment=%s\n',sprintf(' %d',cluster_assignment)));

		fprintf(fp,output_str); 
		fclose(fp);
		
		% In the mat file save the variables which can then be used to
        % create the permutation with the main purpose of being able to plot
        % the points if required.  and reconstruct the information without
        % having to run the algo again. 
		
		
		%quick and dirty hack 
		save(outFileName, 'dataFile','cluster_algo','cluster_assignment','variation_of_information','clustering_error_true','adjRI');

	  end 
	  
	  fprintf(' CE mean=%.4f std=%.4f \n',mean(ce), std(ce)); 
	  fprintf(' VI mean=%.4f std=%.4f \n',mean(vi), std(vi)); 
	  fprintf(' ARI mean=%.4f std=%.4f \n',mean(ari), std(ari)); 
	end
  end
  fprintf('FINISHED AT %.2fs\n',toc); 
  if (plot_points)
	close; 
	plot2Dpoints_with_clusters(points,cluster_assignment);
	title(strcat(dataFile,' | ',cluster_algo,' | K=',sprintf('%d',k),' | \sigma= ',sigma_str(sigma)));
  end
  
