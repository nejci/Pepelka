% adds the Spectral directories to the path. (Ack: modelled on Kevin
% Murphy's BNT add_BNT_to_path)

global SPECTRAL_HOME
global ADDED_PATH

SPECTRAL_HOME=pwd();

%it is global - so we can remove from path all subdirectories
ADDED_PATH=genpath(SPECTRAL_HOME);
addpath(ADDED_PATH);

