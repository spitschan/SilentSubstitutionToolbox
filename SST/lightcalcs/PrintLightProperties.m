function PrintLightProperties(out)


%% Report on stimulus
fprintf('\n');
%fprintf('  * Stimulus is maximum cal file says OneLight can produce\n');
%fprintf('  * Stimulus diameter mm %0.1f, degrees %0.1f\n',2*stimulusRadiusMm,2*stimulusRadiusDeg);
fprintf('  * Stimulus radiance \t%0.1f [log10 watts/[m2-sr]], %0.1f [log10 watts/[cm2-sr]]\n',out.log10SumRadianceWattsPerM2Sr,out.log10SumRadianceWattsPerCm2Sr);
fprintf('  * Stimulus luminance \t%0.1f [candelas/m2\n',out.photopicLuminanceCdM2);
fprintf('  * Stimulus chromaticity x=%0.4f, y=%0.4f\n',out.chromaticityXY(1), out.chromaticityXY(2));
fprintf('  * Stimulus %0.0f (check val %0.0f) scotopic trolands, %0.0f photopic trolands (check val %0.0f)\n',out.irradianceScotTrolands,out.irradianceScotTrolands_check,...
    out.irradiancePhotTrolands,out.irradiancePhotTrolands_check);
fprintf('  * Stimulus %0.1f log10 scotopic trolands, %0.1f log10 photopic trolands\n',log10(out.irradianceScotTrolands),log10(out.irradiancePhotTrolands));
fprintf('  * Stimulus retinal irradiance %0.1f log10 watts/cm2\n',out.log10SumIrradianceWattsPerCm2);
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[cm2-sec]\n',out.log10SumIrradianceQuantaPerCm2Sec);
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[deg2-sec]\n',out.log10SumIrradianceQuantaPerDeg2Sec);
fprintf('  * Stimulus corneal irradiance %0.1f log10 watts/cm2\n',out.log10SumCornealIrradianceWattsPerCm2);
fprintf('  * Stimulus corneal irradiance %0.1f log10 quanta/[cm2-sec]\n',out.log10SumCornealIrradianceQuantaPerCm2Sec);
fprintf('  * Pupil area times LMS: %0.2f, %0.2f, %0.2f\n',...
        theLMSTimesPupilArea(1),theLMSTimesPupilArea(2),theLMSTimesPupilArea(3));