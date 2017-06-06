function [T_cones_sp10,S_cones_sp10] = GetSmithPokornyCones10(S,normalize)
% [T_cones_sp10,S_cones_sp10] = GetSmithPokornyCones10(S,[normalize])
%
% Cao et al. write that they use the spectral sensitivities from Shapiro et
% al. (1996). They in turn use the 1964 10° color-matching functions and then apply
% a transformation to map to cone fundamentals to it. This transformation
% matrix is given in Table 5. Below, this transformation matrix is called M.
%
% Optional argument normalize is false by default.  If set to true, each
% cone sensitivity is normalized to a maximumu of 1.
%
% 10/24/15  dhb  Made it a function.

% Set normalization
if (nargin < 2 | isempty(normalize))
    normalize = false;
end

% Set M matrix reconstruct the photoreceptor sensitivities they are using.
M = [0.15516 0.54308 -0.03287 ; -0.15516 0.45692 0.03287 ; 0 0 1];

% Load the 1964 10° CMFs
wls = (380:1:780)';
load T_xyz1964

% Spline as necessary
if (nargin < 1 | isempty(S))
    S = S_xyz1964;
else
    T_xyz1964 = SplineCmf(S_xyz1964,T_xyz1964,S);
    S_xyz1964 = S;
end

T_cones_sp10 = M*T_xyz1964;
S_cones_sp10 = S;

% Optional normalize
if (normalize)
    for i = 1:size(T_cones_sp10,1)
        T_cones_sp10(i,:) = T_cones_sp10(i,:)/max(T_cones_sp10(i,:));
    end
end

end