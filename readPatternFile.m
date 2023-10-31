function [el, I, t, header, Itotal]=readPatternFile(inputfile)
%READPATTERNFILE  Reads pattern file and returns related info.
%
% Syntax:
%   [el, I, t, header, Itotal] = mats.readPatternFile(inputfile);
%
% Inputs:
%   inputfile - Filename of .pat file.
%   
% Output:
%   el - Electrode identifier/index
%   I  - Current, matched to elements of `el`
%   t  - Time for pulse, matched with elements of `I`
%   header - Metadata parsed from comments at start of file
%   Itotal - Total current, which is half the sum of all rectified current
%
% See also: Contents

arguments
    inputfile (1,1) string % Filename of .pat file
end

f=fopen(inputfile,'r');
header=struct;
el=[];
I=[];
t=[];
elpattern = struct;
while ~feof(f)
    line=fgetl(f);
    if ~isempty(line)
    if line(1)=='#'
        tmp=strsplit(strtrim(strip(line,'#')), ':');
        switch tmp{1}
            case 'Max Density'
                tmp2 = strsplit(strtrim(tmp{2}), ' ');
                header.max_density = struct('value', str2double(tmp2{1}), 'unit', tmp2{2});
            case 'Max Current'
                tmp2 = strsplit(strtrim(tmp{2}), ' ');
                header.max_current = struct('value', str2double(tmp2{1}), 'unit', tmp2{2});
            case 'Location of the Maximum'
                tmp2 = strsplit(tmp{2}, ', ');
                tmp3 = strsplit(tmp2{1}, '[');
                tmp3 = strsplit(tmp3{2}, ' ');
                tmp4 = strsplit(tmp2{2}, ']');
                tmp4 = strsplit(tmp4{1}, ' ');
                header.loc = struct('x', str2double(tmp3{1}), 'y', str2double(tmp4{1}), 'unit', tmp3{2});
            case 'Total Current'
                tmp2 = strsplit(strtrim(tmp{2}), ' ');
                header.total_current = struct('value', str2double(tmp2{1}), 'unit', tmp2{2});
            case 'Area above 80%'
                tmp2 = strsplit(strtrim(tmp{2}), ' ');
                header.area_above_80_pct = struct('value', str2double(tmp2{1}), 'unit', tmp2{2});
            otherwise
                warning("Ignoring unrecognized header entry: <strong>%s</strong>\n", line);
        end
%         header.(strrep(tmp{1}, ' ', '_')) = tmp{2};        
    else

        linef=sscanf(line,'%d\t%g');
        el(end+1)=linef(1);
        I(end+1)=linef(2);
        try
            t(end+1)=linef(3);
        end

    end
    end
end
fclose(f);
if isfield(header, 'total_current')
    Itotal = header.total_current.value;
elseif isfield(header, 'max_current')
    Itotal = header.max_current.value .* sum(I(I > 0));
else
    Itotal = sum(I(I > 0));
end
if isfield(header, 'max_current') && endsWith(inputfile, ".txt")
    I = I .* 1e-5 .* header.max_current.value;
end
el(el > 65) = 65;

end