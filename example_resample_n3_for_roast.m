%EXAMPLE_RESAMPLE_N3_FOR_ROAST  Shows how to resample pre-segmented files for running roast model on ROI.
%
% NOTE: After completing this step, must re-export the patch coordinates
% using the `export_ras_electrodes_2_ijk_custom_locations.m` script in
% NHP-Patch-Experiments repo!

close all force;
clear;
clc;

%% 1. Pick which files to resample.
SUBJ = 'forrest';
NIFTI_TO_RESAMPLE = { ...
    sprintf('example/%s_orig.nii',SUBJ) ...             % This is the T1 (in this case) structural MRI
};
SEGMENTATION_TO_RESAMPLE = { ...
        sprintf('example/%s_T1orT2_masks_orig.nii',SUBJ) ... % This is the segmented version
};

%% 2. Define the ROI (mm) and voxel resolution.
DESIRED_ROI = [ ...
	-26.0,    -25.0,    55.0; ...    % R A S (mm; lower corner)
	 48.0,     10.0,    85.0  ...    % R A S (mm; "upper" corner)
 ];
DESIRED_VOXEL_SIZE = 0.05;           % mm (isotropic voxels)

%% 3. Do the resampling.
% 3a: Resample NIFTI using spm_volume and spm_reslice
for ii = 1:numel(NIFTI_TO_RESAMPLE)
    resample_and_change_bb( ... % Just use spm_volume and spm_reslice
        NIFTI_TO_RESAMPLE{ii},  ...
        DESIRED_VOXEL_SIZE,  ...
        DESIRED_ROI, ...
        "mri", 'Tag', 'resampled' ...
    );
end

% 3b: Resample segmentation using load_nifti_untouch, interp3, and save_untouch_nifti
for ii = 1:numel(SEGMENTATION_TO_RESAMPLE)
    resample_and_change_bb( ... 
        SEGMENTATION_TO_RESAMPLE{ii},  ...
        DESIRED_VOXEL_SIZE,  ...
        DESIRED_ROI, ...
        "segmentation", 'Tag', 'resampled' ...
    );
end

%% 4. Fix the labeling on the output mask.
% % % (Fix file names etc. first) % % %
remap_masks_n3_2_spm(sprintf('example/%s_T1orT2_masks_orig_resampled.nii', SUBJ));
% remap_masks_n3_2_spm(sprintf('example/%s_T1orT2_masks_orig_resampled_remapped.nii',SUBJ), ...
%   'MapInput', struct('gray', 1, 'white', 2, 'csf', 3, 'bone', 4, 'scalp', 5, 'air', 6), ...
%   'MapOutput', struct('gray', 1, 'white', 2, 'csf', 3, 'bone', 4, 'scalp', 5, 'air', 0));

if exist(sprintf('example/%s.nii',SUBJ),'file')~=0
    delete(sprintf('example/%s.nii',SUBJ));
end
if exist(sprintf('example/%s_T1orT2_masks.nii',SUBJ),'file')~=0
    delete(sprintf('example/%s_T1orT2_masks.nii',SUBJ));
end
copyfile(sprintf('example/%s_orig_resampled.nii',SUBJ), sprintf('example/%s.nii', SUBJ), 'f');
copyfile(sprintf('example/%s_T1orT2_masks_orig_resampled_remapped.nii',SUBJ), sprintf('example/%s_T1orT2_masks.nii',SUBJ), 'f');
% copyfile(sprintf('example/%s_T1orT2_masks_orig_resampled_remapped_remapped.nii',SUBJ), sprintf('example/%s_T1orT2_masks.nii',SUBJ), 'f');
delete(sprintf('example/%s_orig_resampled.nii',SUBJ));
delete(sprintf('example/%s_T1orT2_masks_orig_resampled.nii',SUBJ));
delete(sprintf('example/%s_T1orT2_masks_orig_resampled_remapped.nii',SUBJ));
% delete(sprintf('example/%s_T1orT2_masks_orig_resampled_remapped_remapped.nii',SUBJ));
