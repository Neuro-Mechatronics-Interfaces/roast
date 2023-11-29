function [recipe, params] = n3_large_elec_pattern(amp_scale, options)
%N3_LARGE_ELEC_PATTERN  Get N3 currents for large "test" patterns and plug into "recipe" for roast
%
% Syntax:
%   [recipe, params] = n3_large_elec_pattern(amp_scale,'Name',value,...);
%
% Inputs:
%   amp_scale - The total-current amplitude (mA)
%
% Options:
%     options.Conductivities struct = struct('white', 0.126, 'gray', 0.3003, 'csf', 1.7921, 'bone', 0.006, 'skin', 0.330, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
%     options.ArrayElectrodes (1,:) double = 1:64; % Index for any electrodes from array to be used from pattern file.
%     options.ArrayElectrodeSize (1,2) double = [0.5 0.1]; % [radius height] units are mm
%     options.DistantElectrodes (1,:) double = 65; % Index for "distant" return electrode
%     options.DistantElectrodeSize (1,2) double = [9 0.1]; % [radius height] units are mm
%     options.MeshOptions struct = struct('radbound', 4, 'angbound', 10, 'distbound', 0.05, 'reratio', 3, 'maxvol', 5);%
%
% Output:
%   recipe - As used in `roast.m`
%   params - Cell array meant to be used like roast(__,params{:}) to send
%               appropriate electrode size and shape parameters to roast.
%   I_tot - (Scaled) total current used for this pattern, in mA.
%   stim_current - The amount of current on each electrode (for reference)
%
% See also: roast

arguments
    amp_scale (1,1) double;
    % Default `roast` parameters:
    % options.Conductivities = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    % Parameters used in N3 forward modeling (to design experiments):
    options.Conductivities struct = struct('white', 0.126, 'gray', 0.3003, 'csf', 1.7921, 'bone', 0.006, 'skin', 0.330, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    options.ArrayReturnFraction (1,1) double {mustBeInRange(options.ArrayReturnFraction, 0, 1)} = 0.5;
    options.ArrayRingElectrodeSize (1,3) double = [4 8 1]; % [innerRadius outerRadius height] units are mm
    options.ArrayCenterElectrodeSize (1,2) double = [2 1]; % [radius height] units are mm
    options.DistantElectrodeSize (1,2) double = [8 1]; % [radius height] units are mm
    % Default roast (v3) mesh options:
    % options.MeshOptions = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);
    % Custom Mesh Options:
    options.MeshOptions struct = struct('radbound',5,'angbound',30,'distbound',0.2,'reratio',3,'maxvol',4);
    options.Tag {mustBeTextScalar} = "v2";
end

recipe = {'custom1', amp_scale, 'custom2', -amp_scale * options.ArrayReturnFraction, 'custom3', -amp_scale*(1-options.ArrayReturnFraction)};
recipe{end} = -(recipe{2} + recipe{4});

params = cell(1,10);
params{1} = 'elecType';
params{2} = {'disc', 'ring', 'disc'};
params{3} = 'elecSize';
params{4} = {options.ArrayCenterElectrodeSize, options.ArrayRingElectrodeSize, options.DistantElectrodeSize};
params{5} = 'conductivities';
params{6} = options.Conductivities;
params{7} = 'meshOptions';
params{8} = options.MeshOptions;
params{9} = 'simulationTag';
params{10} = sprintf('%s_%d-%d', options.Tag, round(100*options.ArrayReturnFraction), round(100*(1-options.ArrayReturnFraction)));

end