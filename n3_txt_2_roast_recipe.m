function [recipe, params, I_tot, stim_current] = n3_txt_2_roast_recipe(txt_file, amp_scale, options)
%N3_TXT_2_ROAST_RECIPE  Get N3 currents and plug into "recipe" for roast
%
% Syntax:
%   [recipe, params, I_tot, stim_current] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%   [__] = n3_txt_2_roast_recipe(__,'Name',value,...);
%
% Example 1:
%   txt_file = "D:\__Assets__\Data\20220503_theoretical\Patterns\HP20\Jsafety_20_x3000um_y1299um.txt";
%   amp_scale = 0.25; % Used in HP20 experiments in 5-3-22 session
%   [recipe, params] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%
%   This would return recipe and associated parameters for all 64 array
%   electrodes from the pattern .txt file as well as the location and
%   current returned on the "distant" ground electrode.
%
% Example 2:
%   [recipe,params,I_tot,stim_current] = n3_txt_2_roast_recipe( ...
%       txt_file, amp_scale, 'NumArrayElectrodes', 32);
%
%   This would only use the first 32 electrodes specified in the pattern
%   file, along with the location and position of the final electrode
%   parsed from the pattern .txt file.
%
% Inputs:
%   txt_file - Name of .txt file containing stimulating current values
%   amp_scale - The total-current amplitude scaling used for this session
%       -> This was the stimulation parameter that was varied during the
%           actual experiments when tuning MEP thresholds.
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
    txt_file {mustBeTextScalar, mustBeFile}
    amp_scale (1,1) double {mustBeInRange(amp_scale, 0, 1)};
    % Default `roast` parameters:
    % options.Conductivities = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    % Parameters used in N3 forward modeling (to design experiments):
    options.Conductivities struct = struct('white', 0.126, 'gray', 0.3003, 'csf', 1.7921, 'bone', 0.006, 'skin', 0.330, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    options.ArrayElectrodes (1,:) double = 1:64; % Index for any electrodes from array to be used from pattern file.
    options.ArrayElectrodeSize (1,2) double = [0.6 0.2]; % [radius height] units are mm
    options.DistantElectrodes (1,:) double = 65; % Index for "distant" return electrode
    options.DistantElectrodeSize (1,2) double = [5.0 0.2]; % [radius height] units are mm
    % Default roast (v3) mesh options:
    % options.MeshOptions = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);
    % Custom Mesh Options:
    options.MeshOptions struct = struct('radbound',5,'angbound',30,'distbound',0.2,'reratio',3,'maxvol',4);
    options.Tag {mustBeTextScalar} = "";
    options.Version (1,1) double {mustBePositive, mustBeInteger} = 2;
end

[el,I] = readPatternFile(txt_file);

[array_electrodes, array_electrode_indices] = match_subset_indices(el, options.ArrayElectrodes);
array_current = I(array_electrode_indices).*amp_scale.*1e3;
if numel(array_current) == 0
    error("Must have at least one 'array' electrode in montage (otherwise what even is the point of this project).");
end

[distant_electrodes, distant_electrode_indices] = match_subset_indices(el, options.DistantElectrodes);
distant_current = I(distant_electrode_indices).*amp_scale.*1e3;
if numel(distant_current) == 0
    error("Must have at least one 'distant' electrode in montage (can have '0mA' current returned on it though).");
end

stim_current = [array_current, distant_current];
stim_electrodes = [array_electrodes, distant_electrodes];

n_elec_total = numel(stim_electrodes);
if abs(sum(stim_current)) > 0
    i_distant_balancer = numel(array_current) + 1;
    stim_current(i_distant_balancer) = stim_current(i_distant_balancer) - sum(stim_current);
end
recipe = cell(1, 2*n_elec_total);

for ii = 1:n_elec_total
    recipe{1,2*ii-1} = sprintf('custom%d', stim_electrodes(ii));
    recipe{1,2*ii} = stim_current(ii);
end
I_tot = sum(stim_current(stim_current > 0));

params = cell(1,10);
params{1} = 'elecType';
params{2} = repmat({'disc'}, 1, n_elec_total);
params{3} = 'elecSize';
params{4} = [repmat({options.ArrayElectrodeSize},1,numel(array_electrodes)), repmat({options.DistantElectrodeSize},1,numel(distant_electrodes))];
params{5} = 'conductivities';
params{6} = options.Conductivities;
params{7} = 'meshOptions';
params{8} = options.MeshOptions;
params{9} = 'simulationTag';
if strlength(options.Tag) == 0
    [~,pattern_tag,~] = fileparts(txt_file);
    amp_tag = round(amp_scale * 1e2);
%     st = strrep(strrep(string(datetime('now')), ' ', 'T'),':','-');
%     params{10} = sprintf('%s_%d__%s',pattern_tag,amp_tag,st);
    params{10} = sprintf('%s_%dpct_v%d',pattern_tag,amp_tag,options.Version);
else
    params{10} = options.Tag;
end


    function [elements, elementIndices] = match_subset_indices(data_indices, target_indices)
        elements = reshape(target_indices, 1, numel(target_indices));
        elementIndices = nan(size(elements));
        for iElement = 1:numel(elements)
            tmpIndex = find(data_indices==elements(iElement), 1, 'first');
            if ~isempty(tmpIndex)
                elementIndices(iElement) = tmpIndex;
            end
        end
        elements(isnan(elementIndices)) = [];
        elementIndices(isnan(elementIndices)) = [];
    end
end