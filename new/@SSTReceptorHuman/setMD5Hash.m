function obj = setMD5Hash(obj, varargin)

if isempty(obj.MD5Hash)
    tmp = DataHash(obj);
    if strcmp(obj.verbosity, 'high')
       fprintf('* Setting hash: <strong>%s</strong>\n', tmp);
    end
    obj.MD5Hash = tmp;
end