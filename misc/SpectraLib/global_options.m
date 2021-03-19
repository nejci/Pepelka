% This file contains all the global initializations 


global CLUSTER_MULTICUT_NORMALIZE_EV;
CLUSTER_MULTICUT_NORMALIZE_EV=0;


% The various groups of algorithms
clear CA1 CA2 CA3 CA4;
%classical algorithms
CA1{1}='single_linkage';
CA1{2}='ward_linkage';

CA1=CA1';


%single stage spectral 
CA2{1}='njw_ward'; 
CA2{2}='njw_kmeans'; 
CA2{3}='njw_anchor'; 
CA2{4}='mcut_ward'; 
CA2{5}='mcut_kmeans'; 
CA2{6}='mcut_anchor'; 
CA2=CA2';

%double stage spectral 
CA3{1}='njw_njw_ward'; 
CA3{2}='njw_njw_kmeans'; 
CA3{3}='njw_mcut_ward'; 
CA3{4}='njw_mcut_kmeans'; 
CA3{5}='mcut_mcut_ward'; 
CA3{6}='mcut_mcut_kmeans'; 
CA3=CA3';

% recursive spectral ... 
CA4{1}='sm_ncut';
CA4{2}='kvv1_ncut';
CA4{3}='kvv2_ncut';
CA4{4}='sm_cond';
CA4{5}='kvv1_cond';
CA4{6}='kvv2_cond';
CA4{7}='sm_gap'; 
CA4=CA4';

% provided so that run_single_experiment can be used with it easily. 
ALL_ALGOS=[CA1 ;CA2; CA4 ;CA3];


% global constants for some of the algorithm 
global KMEANS_THRESHOLD; 
KMEANS_THRESHOLD=1E-7; 
global KMEANS_MAX_ITER; 
KMEANS_MAX_ITER=100; 
global ANCHOR_NMIN; 
ANCHOR_NMIN=3; 
global SPECTRAL_SIGMA; 
SPECTRAL_SIGMA=0.2; 

% global color vector
% helpful in plotting. 
global DEFAULT_COLOR_VECTOR
DEFAULT_COLOR_VECTOR=['ro'; 'bx'; 'ms'; 'g+'; 'c^'; 'kd'; 'r>'; 'b<' ; 'm*' ; 'gv'; 'cp' ; 'kh'; 
					'bo';'mx';'rs';
					'bs';'mo';'rx';
					'r^';'bv';'md';
					'rv';'bd';'m^';
					'rd';'b^';'mv';
					'go';'gx';'gs';'g^';'gv'; 'gd']; 
DEFAULT_COLOR_VECTOR=strcat(DEFAULT_COLOR_VECTOR,'--'); 

