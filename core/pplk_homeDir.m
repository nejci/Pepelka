function pplkPath = pplk_homeDir(subFolder)
% Returns path to subfolder in Pepelka root directory.
% Default subFolder is core.

%%%FIELD-PATH%%%
pplkPath = 'D:\Nejc\sluzba\LASPP\raziskave\Pepelka\v2';
%%%FIELD-PATH%%%

if ~exist('subFolder','var') || isempty(subFolder)
    pplkPath = [pplkPath,filesep,'core',filesep];
else
    pplkPath =  [pplkPath,filesep,subFolder,filesep];
end