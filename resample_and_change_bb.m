function resample_and_change_bb(mri, desired_voxsize, desired_bb)
%RESAMPLE_AND_CHANGE_BB  Resamples MRI and also changes bounding box (i.e. to restrict ROI to specific part of image)
%
% Syntax:
%   resample_and_change_bb(mri, desired_voxsize, desired_bb);
%
% Inputs:
%     mri {mustBeTextScalar, mustBeFile} - Name of image to resample.
%     desired_voxsize (1,1) double - The desired (isotropic) voxel resolution (mm).
%     desired_bb (2,3) double - The desired bounding box (mm). 
%
% See also: Contents, spm_reslice, spm_vol, resampMesh, roast

arguments
    mri {mustBeTextScalar, mustBeFile};
    desired_voxsize (1,1) double
    desired_bb (2,3) double
end

sz_um = round(diff(desired_bb,1,1).*1e3);
vsz_um = round(desired_voxsize.*1e3);

V = spm_vol(mri);
VV(1:2) = V;
VV(1).mat = spm_matrix([desired_bb(1,:) 0 0 0 ones(1,3).*desired_voxsize])*spm_matrix([-1 -1 -1]);
VV(1).dim = ceil(VV(1).mat \ [desired_bb(2,:) 1]' - 0.1)';
VV(1).dim = VV(1).dim(1:3);
spm_reslice(VV,struct('mean',false,'which',1,'interp',7,'prefix',sprintf('_%dumx%dumx%dum_%dum',sz_um(1),sz_um(2),sz_um(3),vsz_um)));


end