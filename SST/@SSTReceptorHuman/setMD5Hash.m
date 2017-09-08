function obj = setMD5Hash(obj)
% setMD5Hash(ob)
%
% Usage:
%     receptorObj.setMD5Hash(obj);
%
% Description:
%     This method calculates an MD5 checksum on the receptor object and
%     assigns it to the MD5Hash field in the object.
%
% Input:
%     obj - The receptorObj (e.g. from @SSTReceptor or @SSTReceptorHuman)
%
% Output:
%     obj - The receptorObj
%
% Optional key/value pairs:
%     None
%
% See also:
%     @SSTReceptorHuman
%
% 7/25/17    ms       Commented.
% 9/8/17     ms .      Updated header comments
if isempty(obj.MD5Hash)
    tmp = DataHash(obj);
    if strcmp(obj.verbosity, 'high')
       fprintf('* Setting hash: <strong>%s</strong>\n', tmp);
    end
    obj.MD5Hash = tmp;
end