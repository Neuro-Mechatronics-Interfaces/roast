function [elecPara, elecName, injectCurrent] = parse_electrode_parameters(recipe, capType, elecType, elecSize, elecOri)
%PARSE_ELECTRODE_PARAMETERS Parse electrode parameters from input recipe.
%
% Syntax:
%   [elecPara, elecName, injectCurrent] = parse_electrode_parameters(recipe, capType, elecType, elecSize, elecOri);
%
% Inputs:
%   recipe - Main recipe input to `roast` 
%   capType - The EEG system that you want to pick from. Valid options:
%               '1020' | '1010' (default) | '1005' | 'BioSemi' | 'EGI'
%              You can also use customized electrode locations you 
%               defined. Just provide the text file that contains the 
%               electrode coordinates.
%   elecType - The shape of electrode. Valid options:
%               'disc' (default) | 'pad' | 'ring'
%                Note you can specify different shapes to different 
%                   electrodes. In other words, you can place different 
%                   types of electrodes at the same time.
%   elecSize - The size of electrode. 
%               All sizes are in the unit of millimeter (mm). For disc electrodes, sizes follow the format of [radius height], and default size is [6mm 2mm]; for pad electrodes, sizes follow the format of [length width height], and default size is [50mm 30mm 3mm]; for ring electrodes, sizes follow the format of [innerRadius outterRadius height], and default size is [4mm 6mm 2mm].
%             If you're placing only one type of electrode (e.g., either 
%               disc, or pad, or ring), you can use a one-row vector to 
%               customize the size, see Example 7. 
%             If you want to control the size for each electrode 
%               separately (provided you're placing only one type of 
%               electrode), you need to specify the size for each electrode
%               correspondingly in a N-row matrix, where N is the number of
%               electrodes to be placed, see Example 8. 
%             If you're placing more than one type of electrodes and also
%               want to customize the sizes, you need to put the size of 
%               each electrode in a 1-by-N cell (put [] for any electrode 
%               that you want to use the default size), where N is the 
%               number of electrodes to be placed, see Example 9.
%


arguments
    recipe
    capType {mustBeMember(capType, {'1020','1010','1005','biosemi','egi'})}
    elecType
    elecSize
    elecOri
end

if any(~strcmpi(recipe,'leadfield'))
    % check recipe syntax
    if mod(length(recipe),2)~=0
        error('Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair.');
    end

    elecName = (recipe(1:2:end-1))';
    injectCurrent = (cell2mat(recipe(2:2:end)))';
    if abs(sum(injectCurrent))>eps
        error('Electric currents going in and out of the head not balanced. Please make sure they sum to 0.');
    end
    if ~iscellstr(elecType) %#ok<ISCLSTR>
        if ~any(strcmpi(elecType,{'disc','pad','ring'}))
            error('Supported electrodes are: ''disc'', ''pad'' and ''ring''.');
        end
    else
        if length(elecType)~=length(elecName)
            error('You want to place more than 1 type of electrodes, but did not tell ROAST which type for each electrode. Please provide the type for each electrode respectively, as the value for option ''elecType'', in a cell array of length equals to the number of electrodes to be placed.');
        end
        for i=1:length(elecType)
            if ~any(strcmpi(elecType{i},{'disc','pad','ring'}))
                error('Supported electrodes are: ''disc'', ''pad'' and ''ring''.');
            end
        end
    end

    if isempty(elecSize)
        if ~iscellstr(elecType) %#ok<ISCLSTR>
            switch lower(elecType)
                case {'disc'}
                    elecSize = [6 2];
                case {'pad'}
                    elecSize = [50 30 3];
                case {'ring'}
                    elecSize = [4 6 2];
            end
        else
            elecSize = cell(1,length(elecType));
            for i=1:length(elecSize)
                switch lower(elecType{i})
                    case {'disc'}
                        elecSize{i} = [6 2];
                    case {'pad'}
                        elecSize{i} = [50 30 3];
                    case {'ring'}
                        elecSize{i} = [4 6 2];
                end
            end
        end
    else
        if ~iscellstr(elecType) %#ok<ISCLSTR>
            if iscell(elecSize)
                warning('Looks like you''re placing only 1 type of electrodes. ROAST will only use the 1st entry of the cell array of ''elecSize''. If this is not what you want and you meant differect sizes for different electrodes of the same type, just enter ''elecSize'' option as an N-by-2 or N-by-3 matrix, where N is number of electrodes to be placed.');
                elecSize = elecSize{1};
            end
            if any(elecSize(:)<=0)
                error('Please enter non-negative values for electrode size.');
            end
            if size(elecSize,2)~=2 && size(elecSize,2)~=3
                error('Unrecognized electrode sizes. Please specify as [radius height] for disc, [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
            end
            if size(elecSize,1)>1 && size(elecSize,1)~=length(elecName)
                error('You want different sizes for each electrode. Please tell ROAST the size for each electrode respectively, in a N-row matrix, where N is the number of electrodes to be placed.');
            end
            if strcmpi(elecType,'disc') && size(elecSize,2)==3
                error('Redundant size info for Disc electrodes. Please enter as [radius height]');
                %             elecSize = elecSize(:,1:2);
            end
            if any(strcmpi(elecType,{'pad','ring'})) && size(elecSize,2)==2
                error('Insufficient size info for Pad or Ring electrodes. Please specify as [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:,1) < elecSize(:,2))
                error('For Pad electrodes, the width of the pad should not be bigger than its length. Please enter as [length width height]');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:,3) < 3)
                error('For Pad electrodes, the thickness should at least be 3 mm.');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:) > 80)
                warning('You''re placing large pad electrodes (one of its dimensions is bigger than 8 cm). For large pads, the size will not be exact in the model because they will be bent to fit the scalp surface.');
            end
            if strcmpi(elecType,'ring') && any(elecSize(:,1) >= elecSize(:,2))
                error('For Ring electrodes, the inner radius should be smaller than outter radius. Please enter as [innerRadius outterRadius height]');
            end
        else
            if ~iscell(elecSize)
                error('You want to place at least 2 types of electrodes, but only provided size info for 1 type. Please provide complete size info for all types of electrodes in a cell array as the value for option ''elecSize'', or just use defaults by not specifying ''elecSize'' option.');
            end
            if length(elecSize)~=length(elecType)
                error('You want to place more than 1 type of electrodes. Please tell ROAST the size for each electrode respectively, as the value for option ''elecSize'', in a cell array of length equals to the number of electrodes to be placed.');
            end
            for i=1:length(elecSize)
                if isempty(elecSize{i})
                    switch lower(elecType{i})
                        case {'disc'}
                            elecSize{i} = [6 2];
                        case {'pad'}
                            elecSize{i} = [50 30 3];
                        case {'ring'}
                            elecSize{i} = [4 6 2];
                    end
                else
                    if any(elecSize{i}(:)<=0)
                        error('Please enter non-negative values for electrode size.');
                    end
                    if size(elecSize{i},2)~=2 && size(elecSize{i},2)~=3
                        error('Unrecognized electrode sizes. Please specify as [radius height] for disc, [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
                    end
                    if size(elecSize{i},1)>1
                        error('You''re placing more than 1 type of electrodes. Please put size info for each electrode as a 1-row vector in a cell array for option ''elecSize''.');
                    end
                    if strcmpi(elecType{i},'disc') && size(elecSize{i},2)==3
                        error('Redundant size info for Disc electrodes. Please enter as [radius height]');
                        %                     elecSize{i} = elecSize{i}(:,1:2);
                    end
                    if any(strcmpi(elecType{i},{'pad','ring'})) && size(elecSize{i},2)==2
                        error('Insufficient size info for Pad or Ring electrodes. Please specify as [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:,1) < elecSize{i}(:,2))
                        error('For Pad electrodes, the width of the pad should not be bigger than its length. Please enter as [length width height]');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:,3) < 3)
                        error('For Pad electrodes, the thickness should at least be 3 mm.');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:) > 80)
                        warning('You''re placing large pad electrodes (one of its dimensions is bigger than 8 cm). For large pads, the size will not be exact in the model because they will be bent to fit the scalp surface.');
                    end
                    if strcmpi(elecType{i},'ring') && any(elecSize{i}(:,1) >= elecSize{i}(:,2))
                        error('For Ring electrodes, the inner radius should be smaller than outter radius. Please enter as [innerRadius outterRadius height]');
                    end
                end
            end
        end
    end

    if ~iscellstr(elecType) %#ok<ISCLSTR>
        if ~strcmpi(elecType,'pad')
            warning('You''re not placing pad electrodes; customized orientation options will be ignored.');
            elecOri = [];
        else
            if iscell(elecOri)
                allChar = 1;
                for i=1:length(elecOri)
                    if ~ischar(elecOri{i})
                        warning('Looks like you''re only placing pad electrodes. ROAST will only use the 1st entry of the cell array of ''elecOri''. If this is not what you want and you meant differect orientations for different pad electrodes, just enter ''elecOri'' option as an N-by-3 matrix, or as a cell array of length N (put ''lr'', ''ap'', or ''si'' into the cell element), where N is number of pad electrodes to be placed.');
                        elecOri = elecOri{1};
                        allChar = 0;
                        break;
                    end
                end
                if allChar && length(elecOri)~=length(elecName)
                    error('You want different orientations for each pad electrode by using pre-defined keywords in a cell array. Please make sure the cell array has a length equal to the number of pad electrodes.');
                end
            end
            if ~iscell(elecOri)
                if ischar(elecOri)
                    if ~any(strcmpi(elecOri,{'lr','ap','si'}))
                        error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                    end
                else
                    if size(elecOri,2)~=3
                        error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                    end
                    if size(elecOri,1)>1 && size(elecOri,1)~=length(elecName)
                        error('You want different orientations for each pad electrode. Please tell ROAST the orientation for each pad respectively, in a N-by-3 matrix, where N is the number of pads to be placed.');
                    end
                end
            end
        end
    else
        if ~iscell(elecOri)
            elecOri0 = elecOri;
            elecOri = cell(1,length(elecType));
            if ischar(elecOri0)
                if ~any(strcmpi(elecOri0,{'lr','ap','si'}))
                    error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                end
                for i=1:length(elecType)
                    if strcmpi(elecType{i},'pad')
                        elecOri{i} = elecOri0;
                    else
                        elecOri{i} = [];
                    end
                end
            else
                if size(elecOri0,2)~=3
                    error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                end
                numPad = 0;
                for i=1:length(elecType)
                    if strcmpi(elecType{i},'pad')
                        numPad = numPad+1;
                    end
                end
                if size(elecOri0,1)>1
                    if size(elecOri0,1)~=numPad
                        error('You want different orientations for each pad electrode. Please tell ROAST the orientation for each pad respectively, in a N-by-3 matrix, where N is the number of pads to be placed.');
                    end
                else
                    elecOri0 = repmat(elecOri0,numPad,1);
                end
                i0=1;
                for i=1:length(elecType)
                    if strcmpi(elecType{i},'pad')
                        elecOri{i} = elecOri0(i0,:);
                        i0 = i0+1;
                    else
                        elecOri{i} = [];
                    end
                end
            end
        else
            if length(elecOri)~=length(elecType)
                error('You want to place another type of electrodes aside from pad. Please tell ROAST the orienation for each electrode respectively, as the value for option ''elecOri'', in a cell array of length equals to the number of electrodes to be placed (put [] for non-pad electrodes).');
            end
            for i=1:length(elecOri)
                if strcmpi(elecType{i},'pad')
                    if isempty(elecOri{i})
                        elecOri{i} = 'lr';
                    else
                        if ischar(elecOri{i})
                            if ~any(strcmpi(elecOri{i},{'lr','ap','si'}))
                                error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                            end
                        else
                            if size(elecOri{i},2)~=3
                                error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                            end
                            if size(elecOri{i},1)>1
                                error('You''re placing more than 1 type of electrodes. Please put orientation info for each pad electrode as a 1-by-3 vector or one of the three keywords ''lr'', ''ap'', or ''si'' in a cell array for option ''elecOri''.');
                            end
                        end
                    end
                else
                    %                     warning('You''re not placing pad electrodes; customized orientation options will be ignored.');
                    elecOri{i} = [];
                end
            end
        end
    end

else

    fid = fopen('./elec72.loc');
    C = textscan(fid,'%d %f %f %s'); fclose(fid);
    elecName = C{4};
    for i=1:length(elecName)
        elecName{i} = strrep(elecName{i},'.','');
    end
    capType = '1010';
    elecType = 'disc';
    elecSize = [6 2];
    elecOri = [];

end
elecPara = struct(...
    'capType',capType, ...
    'elecType',elecType,...
    'elecSize',elecSize,...
    'elecOri',elecOri);

end