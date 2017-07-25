function obj = setMD5Hash(obj, varargin)
% obj = setMD5Hash(obj, varargin)
%
% This method calculates and sets an MD5 data hash. May be at some point
% useful to store the receptor object as well as a check sum.
%
% 7/25/17   ms      Commented.

if isempty(obj.MD5Hash)
    tmp = DataHash(obj);
    if strcmp(obj.verbosity, 'high')
       fprintf('* Setting hash: <strong>%s</strong>\n', tmp);
    end
    obj.MD5Hash = tmp;
end