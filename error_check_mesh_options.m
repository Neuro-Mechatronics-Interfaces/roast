function meshOptions = error_check_mesh_options(meshOptions, options)
%ERROR_CHECK_MESH_OPTIONS Check input conductivity fields and potentially warn user.
%
% Syntax:
%   meshOptions = error_check_mesh_options(meshOptions, 'Name', value, ...)
%
% Inputs:
%     meshOptions (1,1) struct - Mesh options struct (for cgalmesh.exe, see
%                                iso2mesh for further details).
%
% Options:
%     'Default' (1,1) struct = struct('radbound', 0.5, 'angbound', 30, 'distbound', 0.3, 'reratio', 3, 'maxvol', 10);
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
    meshOptions
    options.Default (1,1) struct = struct('radbound', 0.5, 'angbound', 30, 'distbound', 0.3, 'reratio', 3, 'maxvol', 10);
    options.SuppressWarning (1,1) logical = true;
end

meshOptNam = fieldnames(meshOptions);
if isempty(meshOptNam) || ~all(ismember(meshOptNam,{'radbound';'angbound';'distbound';'reratio';'maxvol'}))
    error('Unrecognized mesh options detected. Supported mesh options are ''radbound'', ''angbound'', ''distbound'', ''reratio'', and ''maxvol''. Please refer to the iso2mesh documentation for more details.');
end
if ~isfield(meshOptions,'radbound')
    meshOptions.radbound = options.Default.radbound;
else
    if ~isnumeric(meshOptions.radbound) || meshOptions.radbound<=0
        error('Please enter a positive number for the mesh option ''radbound''.');
    end
end
if ~isfield(meshOptions,'angbound')
    meshOptions.angbound = options.Default.angbound;
else
    if ~isnumeric(meshOptions.angbound) || meshOptions.angbound<=0
        error('Please enter a positive number for the mesh option ''angbound''.');
    end
end
if ~isfield(meshOptions,'distbound')
    meshOptions.distbound = options.Default.distbound;
else
    if ~isnumeric(meshOptions.distbound) || meshOptions.distbound<=0
        error('Please enter a positive number for the mesh option ''distbound''.');
    end
end
if ~isfield(meshOptions,'reratio')
    meshOptions.reratio = options.Default.reratio;
else
    if ~isnumeric(meshOptions.reratio) || meshOptions.reratio<=0
        error('Please enter a positive number for the mesh option ''reratio''.');
    end
end
if ~isfield(meshOptions,'maxvol')
    meshOptions.maxvol = options.Default.maxvol;
else
    if ~isnumeric(meshOptions.maxvol) || meshOptions.maxvol<=0
        error('Please enter a positive number for the mesh option ''maxvol''.');
    end
end
if any([meshOptions.radbound-options.Default.radbound, ...
        meshOptions.angbound-options.Default.angbound, ...
        meshOptions.distbound-options.Default.distbound, ...
        meshOptions.reratio-options.Default.reratio, ...
        meshOptions.maxvol-options.Default.maxvol]~=0)
    if ~options.SuppressWarning
        warning('You''re changing the advanced meshing options of ROAST. Unless you know what you''re doing, please keep mesh options default.');
    end
end

end