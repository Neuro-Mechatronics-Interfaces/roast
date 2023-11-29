function resample_and_change_bb(mri, desired_voxsize, desired_bb, image_type, options)
%RESAMPLE_AND_CHANGE_BB  Resamples MRI and also changes bounding box (i.e. to restrict ROI to specific part of image)
%
% Syntax:
%   resample_and_change_bb(mri, desired_voxsize, desired_bb, image_type);
%
% Inputs:
%     mri {mustBeTextScalar, mustBeFile} - Name of image to resample.
%     desired_voxsize (1,1) double - The desired (isotropic) voxel resolution (mm).
%     desired_bb (2,3) double - The desired bounding box (mm). 
%     image_type (1,1) string {mustBeTextScalar, mustBeMember(image_type, ["mri", "segmentation"])} = "mri";
%
% See also: Contents, spm_reslice, spm_vol, resampMesh, roast

arguments
    mri {mustBeTextScalar, mustBeFile};
    desired_voxsize (1,1) double
    desired_bb (2,3) double
    image_type (1,1) string {mustBeTextScalar, mustBeMember(image_type, ["mri", "segmentation"])} = "mri";
    options.Tag {mustBeTextScalar} = "";
end

sz_um = round(diff(desired_bb,1,1).*1e3);
vsz_um = round(desired_voxsize.*1e3);
if strcmpi(image_type, "mri")
    V = spm_vol(mri);
    VV(1:2) = V;
    VV(1).mat = spm_matrix([desired_bb(1,:) 0 0 0 ones(1,3).*desired_voxsize])*spm_matrix([-1 -1 -1]);
    VV(1).dim = round(VV(1).mat \ [desired_bb(2,:) 1]')';
    VV(1).dim = VV(1).dim(1:3);
    if strlength(options.Tag) > 0
        pfx = sprintf('_%s', options.Tag);
    else
        pfx = sprintf('_%dumx%dumx%dum_%dum',sz_um(1),sz_um(2),sz_um(3),vsz_um);
    end
    spm_reslice(VV,struct('mean',false,'which',1,'interp',7,'prefix',pfx));
else
    [p,ref_vol_name,e] = fileparts(mri);
    ref_vol_name = strsplit(ref_vol_name,'_');
    ref_vol_name = strcat(ref_vol_name{1},'_orig');
    info = niftiinfo(fullfile(p,strcat(ref_vol_name,e)));
%     info2 = niftiinfo(fullfile(p,sprintf('%s_%dumx%dumx%dum_%dum%s', ref_vol_name, sz_um(1),sz_um(2),sz_um(3),vsz_um, e)));
    q = desired_voxsize ./ info.PixelDimensions;

    M = inv(info.Transform.T');
    ijk_bb = round((M * [desired_bb'; 1,1])'); %#ok<MINV> 
    ijk_bb = min(max(ijk_bb(:,1:3),ones(2,3)),repmat(info.ImageSize,2,1)); % Now we have voxel coordinates, in original space.
    
    vec_I = ijk_bb(1,1):(ijk_bb(2,1));
    vec_J = ijk_bb(1,2):(ijk_bb(2,2));
    vec_K = ijk_bb(1,3):(ijk_bb(2,3)+1);
    [I,J,K] = meshgrid(vec_J, vec_K, vec_I);
    [Iq,Jq,Kq] = meshgrid(ijk_bb(1,2):q(2):ijk_bb(2,2), ...
                ijk_bb(1,3):q(3):(ijk_bb(2,3)+1), ...
                ijk_bb(1,1):q(1):ijk_bb(2,1));
    data = load_untouch_nii(mri);
    img = permute(data.img,[3 2 1]);
    data.img = uint8(interp3(I,J,K,img(vec_K,vec_J,vec_I),Iq,Jq,Kq,'nearest'));
    [p,f,e] = fileparts(mri);
    data.fileprefix = fullfile(p,strrep(f, '_orig', ''));
    data.hdr.dime.dim(2:4) = size(data.img);
    data.hdr.dime.pixdim = [1 ones(1,3).*desired_voxsize 0 0 0 0];
    data.hdr.dime.bitpix = 8;
    data.hdr.hist.descrip = 'spm - realigned';
    data.hdr.hist.qform_code = 2;
    data.hdr.hist.sform_code = 2;
    data.hdr.hist.quatern_b = 0;
    data.hdr.hist.quatern_c = 0;
    data.hdr.hist.quatern_d = 0;
    data.hdr.hist.qoffset_x = desired_bb(1,1);
    data.hdr.hist.qoffset_y = desired_bb(1,2);
    data.hdr.hist.qoffset_z = desired_bb(1,3);
    data.hdr.hist.srow_x = [desired_voxsize 0 0 desired_bb(1,1)];
    data.hdr.hist.srow_y = [0 desired_voxsize 0 desired_bb(1,2)];
    data.hdr.hist.srow_z = [0 0 desired_voxsize desired_bb(1,3)];
    if strlength(options.Tag) > 0
        f = sprintf('%s_%s', f, options.Tag);
    else
        f = sprintf('%s_%dumx%dumx%dum_%dum',f,sz_um(1),sz_um(2),sz_um(3),vsz_um);
    end
    save_untouch_nii(data,fullfile(p,strcat(f,e)));
end

end