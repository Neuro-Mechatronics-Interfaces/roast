function remap_masks_n3_2_spm(mri_masks, options)
%REMAP_MASKS_N3_2_SPM  Fix Slicer3D masks export to have expected TPM indexing
%
% Syntax:
%   remap_masks_n3_2_spm(mri_masks, 'Name', value, ...);
%
% Inputs:
%   mri_masks - The masks.nii file for desired MRI.
%
% See also: Contents
arguments
    mri_masks {mustBeTextScalar, mustBeFile}
    options.MapInput  (1,1) struct = struct('gray', 6, 'white', 2, 'csf', 3, 'bone', 4, 'scalp', 7, 'air', 8);
    options.MapOutput (1,1) struct = struct('gray', 1, 'white', 2, 'csf', 3, 'bone', 4, 'scalp', 5, 'air', 6);
end

data = load_untouch_nii(mri_masks);
MASKS = ["gray", "white", "csf", "bone", "scalp", "air"];
idx = struct;
for ii = 1:6
    idx.(MASKS(ii)) = data.img == options.MapInput.(MASKS(ii));
end
new_image = zeros(size(data.img));
for ii = 1:6
    new_image(idx.(MASKS(ii))) = options.MapOutput.(MASKS(ii));
end
data.img = new_image;

[p,f,e] = fileparts(mri_masks);
f = sprintf('%s_remapped', f);
save_untouch_nii(data,fullfile(p,strcat(f,e)));


end