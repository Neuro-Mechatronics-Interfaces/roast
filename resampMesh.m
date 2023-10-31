function [mriRS,isRS] = resampMesh(mri,isRS,voxsiz)
% [mriRS,isRS] = resampMesh(mri,isRS,voxsiz)
%
% Resample the MRI into `voxsiz` isotropic resolution.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
%   Modified from resampToOneMM by Max Murphy October 2023

vsz = round(voxsiz(1)*1e3)*1e-3;
voxsiz = repmat(vsz,1,3); % Make it isotropic

if isRS
    V = spm_vol(mri);
    [dirname,baseFilename,ext] = fileparts(mri);
    mriRS = [dirname filesep baseFilename sprintf('_%dum', round(vsz*1e3)) ext];
    if exist(mriRS,'file')
        
        warning('%s has already been resampled to %d micron resolution and saved as %s. ROAST will use that file as the input.', mri, round(vsz*1e3), mriRS);
        
    else
        
        fprintf(1,'Resampling %s to %d micron isotropic resolution...', mri, round(vsz*1e3));
        
        bb = spm_get_bbox(V);
        VV(1:2) = V;
        VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
        VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
        VV(1).dim = VV(1).dim(1:3);
        spm_reslice(VV,struct('mean',false,'which',1,'interp',7,'prefix',sprintf('_%dum',round(vsz*1e3))));
        % 'interp' option: 1 for linear, 7 is highest degree (most accurate, slowest)
        % keep using 'prefix' as it's bad to hack SPM variable names
        
        fprintf(1,'%s has been resampled to %d micron isotropic resolution , and is saved as:\n', mri, round(vsz*1e3));
        disp(mriRS);
        disp('It''ll be used as the input for ROAST.');
        
    end
    
else
    
    mriRS = mri;
    
end

end