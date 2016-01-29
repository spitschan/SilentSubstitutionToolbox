function contrasts = ComputeAndReportContrastsFromSpds(string,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,print)
% contrasts = ComputeAndReportContrastsFromSpds(string,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,print)
%
% Report out contrasts function. Assumes that the modulationSpd is the spd
% measured around the background.
%
% See also ComputeAndReportContrastsFromOLPrimaries.

if (nargin < 8 | isempty(print))
    print = true;
end

backgroundReceptors = T_receptors*backgroundSpd;
modulationReceptors = T_receptors*(modulationSpd-backgroundSpd);
contrasts = modulationReceptors ./ backgroundReceptors;
if (print)
    for j = 1:size(T_receptors,1)
        fprintf('\t%s, %s: contrast = %0.1f%%\n',string,photoreceptorClasses{j},100*contrasts(j));
    end
end
