function pool_init(num_threads)
% function for creating and configuring parallel pool
% If pool with desired number of threads (workers) already exists, do
% nothing.
if ~exist('num_threads','var') || isempty(num_threads)
    num_threads = str2double(getenv('NUMBER_OF_PROCESSORS'))-2;
end

createNew = 0;
p = gcp('nocreate');
if isempty(p)
    % no pool here
    createNew = 1;
else
    pSize = p.NumWorkers;
    if pSize ~= num_threads
        createNew = 1;
    end
end

% turn off the local mpiexec implementation
distcomp.feature('LocalUseMpiexec',false); 

if createNew
    delete(gcp('nocreate')); % Delete current parallel pool
    myCluster = parcluster('local'); % load the default profile and edit it
    myCluster.NumWorkers = num_threads; % set number of workers (threads)
    %saveProfile(myCluster); % if you want to save changes
    parpool(myCluster,num_threads);
else
    fprintf(1,'Parpool with desired number of workers already exists.\n');
end


