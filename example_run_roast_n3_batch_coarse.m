%EXAMPLE_RUN_ROAST_N3_BATCH_COARSE  Run `roast` on batch grid for *coarse* (3-elec) N3 pattern
clear;
close all force;
clc;

MRI = 'example/forrestLarge.nii';
GROUP = {'HP10', 'HP20'}; % Can be 'HP10', 'HP20'
ELECTRODE_LOCATIONS_TAG_EXPR = 'customLocations_%d';

for iGroup = 1:2
    % For "large" test model with larger electrode meant to represent the small
    % ones, it would go this way:
    switch GROUP{iGroup}
        case 'HP10'
            Amplitude_mA = 70; % Based on random sweeps @23%-max(~310mA) intensity on Forrest_2022_11_08 
            ArrayReturnFraction = 0.725; % "Jsafety_10_x1000um_y1299um" : 72.4% on array
        case 'HP20'
            Amplitude_mA = 280; % Based on random sweeps @50%-max(~561mA) intensity on Forrest_2022_11_08
            ArrayReturnFraction = 0.95;  % "Jsafety_20_x1000um_y1299um" : 94.19% on array
        otherwise
            error("Unexpected GROUP value: %s", GROUP{iGroup});
    end
    fprintf(1,'Running <strong>roast</strong> model for %s\n', MRI);
    fprintf(1,'\t\t(%4.2f x max-current | %4.1f%% returned on array)\n\n', Amplitude_mA, round(ArrayReturnFraction*100,1));
    
    for iPattern = 1:52
        [recipe, params] = n3_large_elec_pattern(Amplitude_mA, ...
            'ArrayReturnFraction', ArrayReturnFraction, ...
            'Tag', sprintf('v3_%s_p%d', GROUP{iGroup}, iPattern));
        roast(MRI, recipe, params{:}, ...
            'suppressMeshParameterWarning', true, ...
            'suppressConductivityParameterWarning', true, ...
            'voxSize', [0.25 0.25 0.25], ...
            'customElectrodesTag', sprintf(ELECTRODE_LOCATIONS_TAG_EXPR,iPattern), ...
            'visualizeResult', false); % Do that later.
        close all force;
    end
end