%EXAMPLE_RUN_ROAST_N3_SINGLE_COARSE  Run `roast` on single *coarse* (3-elec) N3 pattern
clear;
close all force;

MRI = 'example/forrestLarge.nii';
GROUP = 'HP20'; % Can be 'HP10', 'HP20'
ELECTRODE_LOCATIONS_TAG = 'customLocations';

% For "large" test model with larger electrode meant to represent the small
% ones, it would go this way:
switch GROUP
    case 'HP10'
        Amplitude_mA = 70; % Based on random sweeps @23%-max(~310mA) intensity on Forrest_2022_11_08 
        ArrayReturnFraction = 0.725; % "Jsafety_10_x1000um_y1299um" : 72.4% on array
    case 'HP20'
        Amplitude_mA = 280; % Based on random sweeps @50%-max(~561mA) intensity on Forrest_2022_11_08
        ArrayReturnFraction = 0.95;  % "Jsafety_20_x1000um_y1299um" : 94.19% on array
    otherwise
        error("Unexpected GROUP value: %s", GROUP);
end
[recipe, params] = n3_large_elec_pattern(Amplitude_mA, ...
    'ArrayReturnFraction', ArrayReturnFraction);
clc;
fprintf(1,'Running <strong>roast</strong> model for %s\n', MRI);
fprintf(1,'\t\t(%4.2f x max-current | %4.1f%% returned on array)\n\n', Amplitude_mA, round(ArrayReturnFraction*100,1));
pause(2.0);
roast(MRI, recipe, params{:}, ...
    'suppressMeshParameterWarning', true, ...
    'suppressConductivityParameterWarning', true, ...
    'voxSize', [0.25 0.25 0.25], ...
    'customElectrodesTag', ELECTRODE_LOCATIONS_TAG);