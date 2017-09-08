function CheckSSTInstalled
% CheckSSTInstalled.m
%
% Check that SST and its dependencies are set up properly.
%
% If ToolboxToolbox (TbTb, http://github.com/ToolboxHub/ToolboxToolbox) is
% installed, this will go ahead and do the required setup for you.  Note
% that the TbTb call does NOT reset the path before it does its stuff, so
% that calling this routine will not clobber a larger setup that is already
% in place.
%
% If you don't have TbTb then this routine prints out an error message
% with info about what you need to install.

% 07/24/17   ms      Written.
% 09/08/17   dhb     Edit so it doesn't reinstall if things are already there.
%            dhb     Edit so that tbUse call doesn't reset Matlab path.
%            dhb     Make this a function, not a script.

% Check if the SilentSubstitutionToolbox is on the path.
% If not and if you have TbTb, then set up.
if ~any(strfind(path, 'SilentSubstitutionToolbox'))
    if any(strfind(path, 'ToolboxToolbox'))
        % Set up paths
        tbUse('SilentSubstitutionToolbox','reset','as-is');
    else
        error('* The SilentSubstitutionToolbox and its dependencies required to run this demo.\nPlease install from https://github.com/spitschan/SilentSubstitutionToolbox and add it and its dependencies to your MATLAB path.');
    end
end

% Check if the BrainardLabToolbox is on the path
if ~any(strfind(path, 'BrainardLabToolbox'))
    if any(strfind(path, 'ToolboxToolbox'))
        tbUse('SilentSubstitutionToolbox','reset','as-is');
    else
        error('* The BrainardLabToolbox is required to run SilentSubstitutionToolbox.\nPlease install from https://github.com/DavidBrainard/BrainardLabToolbox and add to your MATLAB path.');
    end
end

% Check if PTB-3 is on the path
if ~any(strfind(path, 'Psychtoolbox-3'))
    if any(strfind(path, 'ToolboxToolbox'))
        tbUse('SilentSubstitutionToolbox','reset','as-is');
    else
        error('Psychtoolbox-3 is required to run SilentSubstitutionToolbox.\nPlease install from https://github.com/Psychtoolbox-3/Psychtoolbox-3 and add to your MATLAB path.');
    end
end

end

