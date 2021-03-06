function sstValidateFullAll(varargin)
% Run unit tests for SST
%
% Syntax:
%    sstValidateFullAll
%
% Description:
%    Run the set of unit tests for the SilentSubstitutionToolbox. Return
%    true or false depending on success. Designed for use with Jenkins
%    build integration.
%

%% Close all figures so that we start with a clean slate
close all;

%% Set return flag
success = true;

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'SilentSubstitutionToolbox';

%% Here are some tests
status = ReceptorIsolateDemo('validate','basichuman','whichDirectionNumber',1);
if (~status)
    fprintf('ReceptorIsolateDemo failed\n');
    success = false;
end

status = SSTReceptorCheckAgreement();
if (~status)
    fprintf('SSTReceptorCheckAgreement failed\n');
    success = false;
end

status = SSTReceptorDemo('validate','basichuman');
if (~status)
    fprintf('SSTReceptorDemo failed\n');
    success = false;
end

%% Close up all figs
close all;

%% Report whether we are OK
assert(success, 'One or more validations failed.');

end
