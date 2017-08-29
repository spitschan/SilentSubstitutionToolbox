function PrintLightProperties(out)
% PrintLightProperties(out)
%
% Print out the properties of the stimulus.


%% Report on stimulus
fprintf('* <strong>Radiance</strong>\n');
fprintf('\t%0.1f log10 watts/[m2-sr]\n',out.log10SumRadianceWattsPerM2Sr);
fprintf('\t%0.1f log10 watts/[cm2-sr]\n', out.log10SumRadianceWattsPerCm2Sr);
fprintf('\n');
fprintf('* <strong>Luminance</strong>\n');
fprintf('\t%0.1f cd/m2\n',out.photopicLuminanceCdM2);
fprintf('\t%0.1f log10 cd/m2\n',log10(out.photopicLuminanceCdM2));
fprintf('\n');
fprintf('* <strong>Illuminance</strong>\n');
fprintf('\t%0.1f sc Td\n', out.irradianceScotTrolands);
fprintf('\t%0.1f log10 sc Td\n', log10(out.irradianceScotTrolands));
fprintf('\t%0.1f ph Td\n', out.irradiancePhotTrolands);
fprintf('\t%0.1f log10 ph Td\n', log10(out.irradiancePhotTrolands));
fprintf('\n');
fprintf('* <strong>Chromaticity</strong>\n');
fprintf('\tx=%0.4f\n\ty=%0.4f\n',out.chromaticityXY(1), out.chromaticityXY(2));
fprintf('\n');
fprintf('* <strong>Retinal irradiance</strong>\n');
fprintf('\t%0.1f log10 watts/cm2\n',out.log10SumIrradianceWattsPerCm2);
fprintf('\t%0.1f log10 quanta/[cm2-sec]\n',out.log10SumIrradianceQuantaPerCm2Sec);
fprintf('\t%0.1f log10 quanta/[deg2-sec]\n',out.log10SumIrradianceQuantaPerDeg2Sec);
fprintf('\n');
fprintf('* <strong>Corneal irradiance</strong>\n');
fprintf('\t%0.1f log10 watts/cm2\n',out.log10SumCornealIrradianceWattsPerCm2);
fprintf('\t%0.1f log10 quanta/[cm2-sec]\n',out.log10SumCornealIrradianceQuantaPerCm2Sec);
fprintf('\n');
fprintf('* <strong>Melanopic irradiance</strong>\n');
fprintf('\t%0.1f log10 quanta/[cm2-sec] <strong>[corneal]</strong>\n',out.log10SumMelCornealIrradianceQuantaPerCm2Sec);
fprintf('\t%0.1f log10 quanta/[cm2-sec] <strong>[retinal]</strong>\n',out.log10SumMelIrradianceQuantaPerCm2Sec);