function setup
% Pepelka - data clustering package
% Setup script:
%   - define absolute path to Pepelka in pplk_homeDir
%   - add Pepelka/core to Matlab path
%   - compile mex files
% -------------------------------------------------------------------------
sep = repmat('-',1,50);
fprintf('%s\nPepelka - data clustering package - SETUP\n%s\n\n',sep,sep);

% Set home folder
fprintf('Setting home folder in "core/pplk_homeDir.m" to:\n"%s" ',pwd());
fid = fopen(['core',filesep,'pplk_homeDir.m'],'r');
s = fread(fid,'*char')';
fclose(fid);
tokens = strsplit(s,'%%%FIELD-PATH%%%\n');
tokens{2} = ['pplkPath = ''',pwd(),''';'];
fid = fopen(['core',filesep,'pplk_homeDir.m'],'w');
fprintf(fid,'%s%%%%%%FIELD-PATH%%%%%%\n%s\n%%%%%%FIELD-PATH%%%%%%\n%s',tokens{1},tokens{2},tokens{3});
fclose(fid);
fprintf('[OK]\n%s\n',sep);

answ = input('Add "core" folder to the Matlab path? ([y]/n) > ','s');
% Add to path
if isempty(answ) || lower(answ)=='y' || strcmpi(answ,'yes')
    addpath('core');
    savepath();
    fprintf('[OK]\n%s\n',sep);
end



% Compile
answ = input('Do you want to compile MEX files required by some methods now? (You should set-up a compiler first.) ([y]/n) > ','s');
if isempty(answ) || lower(answ)=='y' || strcmpi(answ,'yes')
    
    fprintf(1,'Compiling MEX files\n');
    
    % DICLENS
    % Note: mst_mex.mex* and components_mex.mex* are from MATLAB BGL package:
    % https://dgleich.github.io/matlab-bgl/ (misc\matlab_BGL\private)
    % 'fast' version
    folder = ['methods',filesep,'DICLENS',filesep,'fast',filesep];
    fprintf(1,'\n%s\n%s\n',sep,folder);
    mex('-largeArrayDims', '-outdir', folder, [folder,'computeP_mex.c']);
    mex('-largeArrayDims', '-outdir', folder, [folder,'computeG_mex.c'], 'OPTIMFLAGS="/openmp $OPTIMFLAGS"');
    mex('-largeArrayDims', '-outdir', folder, [folder,'DiclensMainLoop_mex.c']);
    mex('-largeArrayDims', '-outdir', folder, [folder,'majorityVoting_mex.c']);
    folder_matlabBGL = ['misc',filesep,'matlab_BGL',filesep,'private',filesep];
    copyfile([folder_matlabBGL,'mst_mex.',mexext],folder);
    copyfile([folder_matlabBGL,'components_mex.',mexext],folder);
    % 'fastModW' version
    folder = ['methods',filesep,'DICLENS',filesep,'fastModW',filesep];
    fprintf(1,'\n%s\n%s\n',sep,folder);
    mex('-largeArrayDims', '-outdir', folder, [folder,'computeP_mex.c']);
    mex('-largeArrayDims', '-outdir', folder, [folder,'computeG_mex.c'], 'OPTIMFLAGS="/openmp $OPTIMFLAGS"');
    mex('-largeArrayDims', '-outdir', folder, [folder,'DiclensMainLoop_mex.c']);
    mex('-largeArrayDims', '-outdir', folder, [folder,'majorityVoting_mex.c']);
    folder_matlabBGL = ['misc',filesep,'matlab_BGL',filesep,'private',filesep];
    copyfile([folder_matlabBGL,'mst_mex.',mexext],folder);
    copyfile([folder_matlabBGL,'components_mex.',mexext],folder);
    
    % LCE
    folderDst = ['methods',filesep,'LCE',filesep];
    folder = [folderDst,'MEX',filesep];
    fprintf(1,'\n%s\n%s\n',sep,folder);
    mex('-largeArrayDims', '-outdir', folderDst, [folder,'cts_S_mex.c'], 'OPTIMFLAGS="/openmp $OPTIMFLAGS"');
    mex('-largeArrayDims', '-outdir', folderDst, [folder,'weightCl.c'], 'OPTIMFLAGS="/openmp $OPTIMFLAGS"');
    
    % SPECLS
    folder = ['methods',filesep,'SPECLS',filesep];
    fprintf(1,'\n%s\n%s\n',sep,folder);
    mexFiles = dir([folder,'*.cpp']);
    mexFiles = {mexFiles.name};
    for i=1:length(mexFiles)
        mex('-outdir', folder, [folder,mexFiles{i}]);
    end
    
    % DNs index
    folderDst = ['validation',filesep,'indexDNmod',filesep];
    folder = [folderDst,'MEX',filesep];
    fprintf(1,'\n%s\n%s\n',sep,folder);
    mex('-largeArrayDims', '-outdir', folderDst, [folder,'gabriel.c']);
        
    fprintf('\n%s\nCompilation [OK]\n%s\n',sep,sep);
end

fprintf('\n%s\nEverything is prepared. Enjoy using Pepelka!\n',sep);
fprintf('Report bugs to nejc.ilc@fri.uni-lj.si. Thanks!\n');
fprintf('To view a demo run pplk_demo.\n');
