function isolatingPrimary = ReceptorIsolateWrapper(mode, T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
% isolatingPrimary = ReceptorIsolateWrapper(mode, T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
%
% Wrapper functionc all different versions of ReceptorIsolate:
%   - ReceptorIsolate(...)   - Standard method
%   - ReceptorIsolateStrictPhotonIncrease(...)  - Only allows for
%                                                 increments in the primaries
%   - ReceptorIsoalteStrictPhotonDecrease(...)  - Only allows for
%                                                 decrements in the primaries
%
% 6/17/2013     spitschan       Written. This should probably be integrated
%                               into the main ReceptorIsolate code somehow.

%% Execute some consistency checks
% Check whether the desired contrasts vector was passed, and if so check consistency of its dimensions.
if (nargin < 11)
    desiredConstrast = [];
end
if ~isempty(desiredContrasts)
    if length(whichReceptorsToIsolate) ~= length(desiredContrasts)
        error('Size of whichReceptorsToIsolate and of desired contrasts vector do not line up')
    end
end

%% Default for ambientSpd
if (nargin < 12 || isempty(ambientSpd))
    ambientSpd = zeros(size(B_primary,1),1);
end

%% Select the call to the respective functions
switch mode
    case 'Standard'
        isolatingPrimary = ReceptorIsolate(T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
    case 'StrictPhotonIncrease'
        isolatingPrimary = ReceptorIsolateStrictPhotonIncrease(T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
    case 'StrictPhotonDecrease'
        isolatingPrimary = ReceptorIsolateStrictPhotonDecrease(T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
    case 'EnforceSpectralChange'
        isolatingPrimary = ReceptorIsolateEnforceSpectralChange(T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrasts,ambientSpd);
    case 'None'
        isolatingPrimary = backgroundPrimary;
    otherwise
        error(['ReceptorIsolate mode "' mode '"not defined']);
end