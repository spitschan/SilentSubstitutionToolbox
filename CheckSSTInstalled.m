% CheckSSTInstalled.m
%
% Check that SST and its dependencies are set up properly.
%
% 7/24/17   ms      Written.

% Check if the TbTb is on the path. If yes, we use tbUse to set up
% dependencies. If not, we alert the user to download
if any(strfind(path, 'ToolboxToolbox'))
    % Set up paths
    r = tbUse('SilentSubstitutionToolbox');
else
    
    % Check if the BrainardLabToolbox is on the path
    if ~any(strfind(path, 'BrainardLabToolbox'))
        error('* The BrainardLabToolbox is required to run SilentSubstitutionToolbox.\nPlease install from https://github.com/DavidBrainard/BrainardLabToolbox and add to your search path.')
    end
    
    % Check if PTB-3 is on the path
    if ~any(strfind(path, 'Psychtoolbox-3'))
        error('Psychtoolbox-3 is required to run SilentSubstitutionToolbox.\nPlease install from https://github.com/Psychtoolbox-3/Psychtoolbox-3 and add to your search path path.')
    end
end