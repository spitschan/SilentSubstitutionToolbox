function outVal = SSTGetDefaultParameterRange(whichParameter, outputType)
% SSTGetDefaultParameterRange(whichParameter)
%
% Usage:
%     rangeVals = SSTGetDefaultParameterRange(whichParameter)
%
% Description:
%     This returns the default range of parameters for the Asano model.
%
% Input:
%     whichParameter - String containing which parameter is requested.
%                      Possible options are:
%                           'dlens' (Default)
%                           'dmac'
%                           'dphotopigmentL'
%                           'dphotopigmentM'
%                           'dphotopigmentS'
%                           'lambdaMaxShiftL'
%                           'lambdaMaxShiftM'
%                           'lambdaMaxShiftS'
%                           'obsPupilDiameterMm'
%     outputType - String determing the ouput. Options are:
%                           'start' - start of range
%                           'end' - end of range
%                           'range' - both start and end
%
% Output:
%     outVal - requested range value (start, end, or full range)
%
% Optional key/value pairs:
%     None.

% 9/8/17  ms  Added header comments.

% Define the parameter ranges
switch whichParameter
    case 'dlens'
        rangeVals = [-50 50];
    case 'dmac'
        rangeVals = [-100 100];
    case 'dphotopigmentL'
        rangeVals = [-50 50];
    case 'dphotopigmentM'
        rangeVals = [-50 50];
    case 'dphotopigmentS'
        rangeVals = [-50 50];
    case 'lambdaMaxShiftL'
        rangeVals = [-5 5];
    case 'lambdaMaxShiftM'
        rangeVals = [-5 5];
    case 'lambdaMaxShiftS'
        rangeVals = [-5 5];
    case 'obsPupilDiameterMm'
        rangeVals = [-2 9];
    otherwise
        error('Unknown parameter type');
end

% Return what is requested
switch outputType
    case 'start'
        outVal = rangeVals(1);
    case 'end'
        outVal = rangeVals(end);
    case 'range'
        outVal = rangeVals;
    otherwise
        error('Unknown range type');
end