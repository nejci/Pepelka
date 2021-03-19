% Compile script
mex -largeArrayDims computeP_mex.c;
mex -largeArrayDims computeG_mex.c OPTIMFLAGS="/openmp $OPTIMFLAGS";
mex -largeArrayDims DiclensMainLoop_mex.c;
mex -largeArrayDims majorityVoting_mex.c;