function pplkPath = pplk_homeDir(subFolder)
% Returns path to subfolder in Pepelka root directory.
% Default subFolder is core.

P = load('pplk_userprefs.mat');
pplkPath = P.homeDir;

if ~exist('subFolder','var') || isempty(subFolder)
    pplkPath = [pplkPath,filesep,'core',filesep];
else
    pplkPath =  [pplkPath,filesep,subFolder,filesep];
end