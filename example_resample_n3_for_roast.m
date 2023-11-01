%EXAMPLE_RESAMPLE_N3_FOR_ROAST  Shows how to resample pre-segmented files for running roast model on ROI.

close all force;
clear;
clc;

%% 1. Pick which files to resample.
NIFTI_TO_RESAMPLE = { ...
    'example/forrest_orig.nii', ...             % This is the T1 (in this case) structural MRI
    'example/forrest_T1orT2_masks_orig.nii' ... % This is the segmented version
};

%% 2. Define the ROI (mm) and voxel resolution.
DESIRED_ROI = [ ...
	-23.0,    -17.0,    57.0; ...    % R A S (mm; lower corner)
	 45.0,      3.0,    74.0  ...    % R A S (mm; "upper" corner)
 ];
DESIRED_VOXEL_SIZE = 0.15;           % mm (isotropic voxels)

%% 3. Do the resampling.
for ii = 1:numel(NIFTI_TO_RESAMPLE)
    resample_and_change_bb( ... % Just use spm_volume and spm_reslice
        NIFTI_TO_RESAMPLE{ii},  ...
        DESIRED_VOXEL_SIZE,  ...
        DESIRED_ROI ...
    );
end