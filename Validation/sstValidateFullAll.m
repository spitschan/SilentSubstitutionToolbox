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

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'SilentSubstitutionToolbox';

%% Put in the tests here

%% Set return flag
success = true;
assert(success, 'One or more validations failed.');

end