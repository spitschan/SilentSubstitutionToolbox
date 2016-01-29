function contrasts = ComputeAndReportContrasts(string,photoreceptorClasses,T_receptors,B_primary,backgroundPrimary,theModulationPrimary,ambientSpd,print)
% contrasts = ComputeAndReportContrasts(string,photoreceptorClasses,T_receptors,B_primary,backgroundPrimary,theModulationPrimary,ambientSpd,print)
% Report out contrasts function

if (nargin < 8 | isempty(print))
    print = true;
end

backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*theModulationPrimary;
contrasts = modulationReceptors ./ backgroundReceptors;
if (print)
    for j = 1:size(T_receptors,1)
        fprintf('\t%s, %s: contrast = %0.1f%%\n',string,photoreceptorClasses{j},100*contrasts(j));
    end
end
