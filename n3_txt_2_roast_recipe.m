function [recipe, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale)
%N3_TXT_2_ROAST_RECIPE  Get N3 currents and plug into "recipe" for roast
%
% Syntax:
%   [recipe, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%
% Example:
%   txt_file = "D:\__Assets__\Data\20220503_theoretical\Patterns\HP20\Jsafety_20_x3000um_y1299um.txt";
%   amp_scale = 0.25; % Used in HP20 experiments in 5-3-22 session
%   [recipe, I_tot] = n3_txt_2_roast_recipe(txt_file, amp_scale);
%
% Inputs:
%   txt_file - Name of .txt file containing stimulating current values
%   amp_scale - The total-current amplitude scaling used for this session
%       -> This was the stimulation parameter that was varied during the
%           actual experiments when tuning MEP thresholds.
%
% Output:
%   recipe - As used in `roast.m`
%   I_tot - (Scaled) total current used for this pattern, in mA.
%
% See also: roast

arguments
    txt_file {mustBeTextScalar, mustBeFile}
    amp_scale (1,1) double {mustBeInRange(amp_scale, 0, 1)};
end

stim_current = zeros(65,1);
[el,I] = readPatternFile(txt_file);
stim_current(el) = I .* amp_scale .* 1e3; % Convert to mA
recipe = cell(1, 2*numel(stim_current));
for ii = 1:numel(stim_current)
    recipe{1,2*ii-1} = sprintf('custom%d', ii);
    recipe{1,2*ii} = stim_current(ii);
end
I_tot = sum(stim_current(stim_current > 0));

end