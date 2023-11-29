function conductivities = error_check_conductivities(conductivities, n_elec, options)
%ERROR_CHECK_CONDUCTIVITIES Check input conductivity fields and potentially warn user.
%
% Syntax:
%   conductivities = error_check_conductivities(conductivities, n_elec);
%   __ = error_check_conductivities(__,'Name', value, ...);
%
% Inputs:
%     conductivities (1,1) struct - conductivities struct
%     n_elec (1,1) double {mustBePositive, mustBeInteger} - # of electrode names
%
% Options:
%     'DefaultConductivity' = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
%     'SuppressWarning' (1,1) logical = true;
%
% Output:
%   conductivities - Same as input struct but checked for errors and
%                    potentially 'gel' and 'electrode' field values have 
%                    been replicated to match total number of electrodes 
%                    (if they were given as scalar).
%
% See also: Contents, roast

arguments
    conductivities (1,1) struct
    n_elec (1,1) double {mustBePositive, mustBeInteger}
    options.DefaultConductivity = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    options.SuppressWarning (1,1) logical = true;
end

conductivitiesNam = fieldnames(conductivities);
if isempty(conductivitiesNam) || ~all(ismember(conductivitiesNam,{'white';'gray';'csf';'bone';'skin';'air';'gel';'electrode'}))
    error('Unrecognized tissue names detected. Supported tissue names in the conductivity option are ''white'', ''gray'', ''csf'', ''bone'', ''skin'', ''air'', ''gel'' and ''electrode''.');
end
if ~isfield(conductivities,'white')
    conductivities.white = options.DefaultConductivity.white;
else
    if ~isnumeric(conductivities.white) || any(conductivities.white(:)<=0)
        error('Please enter a positive number for the white matter conductivity.');
    end
    if length(conductivities.white(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'gray')
    conductivities.gray = options.DefaultConductivity.gray;
else
    if ~isnumeric(conductivities.gray) || any(conductivities.gray(:)<=0)
        error('Please enter a positive number for the gray matter conductivity.');
    end
    if length(conductivities.gray(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'csf')
    conductivities.csf = options.DefaultConductivity.csf;
else
    if ~isnumeric(conductivities.csf) || any(conductivities.csf(:)<=0)
        error('Please enter a positive number for the CSF conductivity.');
    end
    if length(conductivities.csf(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'bone')
    conductivities.bone = options.DefaultConductivity.bone;
else
    if ~isnumeric(conductivities.bone) || any(conductivities.bone(:)<=0)
        error('Please enter a positive number for the bone conductivity.');
    end
    if length(conductivities.bone(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'skin')
    conductivities.skin = options.DefaultConductivity.skin;
else
    if ~isnumeric(conductivities.skin) || any(conductivities.skin(:)<=0)
        error('Please enter a positive number for the skin conductivity.');
    end
    if length(conductivities.skin(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'air')
    conductivities.air = options.DefaultConductivity.air;
else
    if ~isnumeric(conductivities.air) || any(conductivities.air(:)<=0)
        error('Please enter a positive number for the air conductivity.');
    end
    if length(conductivities.air(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
end
if ~isfield(conductivities,'gel')
    conductivities.gel = options.DefaultConductivity.gel;
else
    if ~isnumeric(conductivities.gel) || any(conductivities.gel(:)<=0)
        error('Please enter a positive number for the gel conductivity.');
    end
    if length(conductivities.gel(:))>1 && length(conductivities.gel(:))~=n_elec
       error('You want to assign different conductivities to the conducting media under different electrodes, but didn''t tell ROAST clearly which conductivity each electrode should use. Please follow the order of electrodes you put in ''recipe'' to give each of them the corresponding conductivity in a vector as the value for the ''gel'' field in option ''conductivities''.');
    end
end
if ~isfield(conductivities,'electrode')
    conductivities.electrode = options.DefaultConductivity.electrode;
else
    if ~isnumeric(conductivities.electrode) || any(conductivities.electrode(:)<=0)
        error('Please enter a positive number for the electrode conductivity.');
    end
    if length(conductivities.electrode(:))>1 && length(conductivities.electrode(:))~=n_elec
       error('You want to assign different conductivities to different electrodes, but didn''t tell ROAST clearly which conductivity each electrode should use. Please follow the order of electrodes you put in ''recipe'' to give each of them the corresponding conductivity in a vector as the value for the ''electrode'' field in option ''conductivities''.');
    end
end
if ~options.SuppressWarning
    if any([conductivities.white-options.DefaultConductivity.white, ...
            conductivities.gray-options.DefaultConductivity.gray, ...
            conductivities.csf-options.DefaultConductivity.csf, ...
            conductivities.bone-options.DefaultConductivity.bone, ...
            conductivities.skin-options.DefaultConductivity.skin, ...
            conductivities.air-options.DefaultConductivity.air, ...
            conductivities.gel-options.DefaultConductivity.gel, ...
            conductivities.electrode-options.DefaultConductivity.electrode]~=0)
        warning('You''re changing the advanced conductivity options of ROAST. Unless you know what you''re doing, please keep conductivity values default.');
    end
end

if length(conductivities.gel(:))==1
    conductivities.gel = repmat(conductivities.gel,1,n_elec);
end
if length(conductivities.electrode(:))==1
    conductivities.electrode = repmat(conductivities.electrode,1,n_elec);
end

end