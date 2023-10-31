function [recipe, params, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale, options)
%N3_TXT_2_ROAST_RECIPE  Get N3 currents and plug into "recipe" for roast
%
% Syntax:
%   [recipe, params, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%
% Example:
%   txt_file = "D:\__Assets__\Data\20220503_theoretical\Patterns\HP20\Jsafety_20_x3000um_y1299um.txt";
%   amp_scale = 0.25; % Used in HP20 experiments in 5-3-22 session
%   [recipe, params, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%
% Inputs:
%   txt_file - Name of .txt file containing stimulating current values
%   amp_scale - The total-current amplitude scaling used for this session
%       -> This was the stimulation parameter that was varied during the
%           actual experiments when tuning MEP thresholds.
%
% Output:
%   recipe - As used in `roast.m`
%   params - Cell array meant to be used like roast(__,params{:}) to send
%               appropriate electrode size and shape parameters to roast.
%   I_tot - (Scaled) total current used for this pattern, in mA.
%
% See also: roast

arguments
    txt_file {mustBeTextScalar, mustBeFile}
    amp_scale (1,1) double {mustBeInRange(amp_scale, 0, 1)};
    options.NumArrayElectrodes (1,1) double = 64;
    options.ArrayElectrodeSize (1,:) double = []; % [innerRadius outerRadius height] units are mm
    options.DistantElectrodeSize (1,:) double = []; % [radius height] units are mm
end

n_elec = options.NumArrayElectrodes + 1; % Add 1 for "distant return" electrode.
stim_current = zeros(n_elec,1);
[el,I] = readPatternFile(txt_file);
stim_current(el) = I .* amp_scale .* 1e3; % Convert to mA
recipe = cell(1, 2*numel(stim_current));
for ii = 1:n_elec
    recipe{1,2*ii-1} = sprintf('custom%d', ii);
    recipe{1,2*ii} = stim_current(ii);
end
I_tot = sum(stim_current(stim_current > 0));

params = cell(1,5);
params{1} = 'elecType';
params{2} = [repmat({'ring'}, 1, options.NumArrayElectrodes), {'disc'}];
params{3} = 'elecSize';
params{4} = [repmat({options.ArrayElectrodeSize},1,options.NumArrayElectrodes), {options.DistantElectrodeSize}];
params{5} = 'conductivities';
% Default `roast` parameters:
% params{6} = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
% Parameters used in N3 forward modeling (to design experiments):
params{6} = struct('white', 0.126, 'gray', 0.3003, 'csf', 1.7921, 'bone', 0.006, 'skin', 0.330, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
params{7} = 'simulationTag';
[~,pattern_tag,~] = fileparts(txt_file);
amp_tag = round(amp_scale * 1e2);
params{8} = sprintf('%s__%d',pattern_tag,amp_tag);

end