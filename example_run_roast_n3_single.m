%EXAMPLE_RUN_ROAST_N3_SINGLE  Run `roast` on single N3 pattern
clear;
close all force;

%% 1. Set parameters.
MRI = 'example/forrest.nii';
PATTERN = "Jsafety_20_x1000um_y1299um";
AMPLITUDE = 0.25;
PATTERN_ROOT = "R:/NMLShare/users/mdmurphy/Projects/N3/SharpFocus/20220503_theoretical/Patterns";
PATTERN_GROUP = "HP20";

%% 2. Run `roast` model.
txt_file = sprintf("%s/%s/%s.txt", PATTERN_ROOT, PATTERN_GROUP, PATTERN);
[recipe, params] = n3_txt_2_roast_recipe(txt_file, AMPLITUDE); 
clc;
fprintf(1,'Running <strong>roast</strong> model for %s using pattern <strong>%s</strong>\n', MRI, PATTERN);
fprintf(1,'\t\t(%4.2f x max-current)\n\n', AMPLITUDE);
pause(2.0);
roast(MRI, recipe, params{:});