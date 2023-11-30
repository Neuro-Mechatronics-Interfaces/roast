function roast(subj,recipe,options)
%ROAST - roast(subj,recipe) - Main function of ROAST.
%
% Syntax:
%   ROAST(subj,recipe,'Name',value,...)
%
% Inputs:
%     subj {mustBeTextScalar} = 'example/MNI152_T1_1mm.nii';
%     recipe (1,:) cell = {'Fp1',1,'P4',-1};
%     options.capType {mustBeMember(options.capType, {'1020','1010','1005','biosemi','egi'})} = '1010';
%     options.elecType = [];
%     options.elecSize = [];
%     options.elecOri = 'lr'; % Can be 'lr', 'ap', 'si', or (1,3) vector giving precise vector of long axis
%     options.T2 = [];
%     options.meshOptions (1,1) struct = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);
%     options.simulationTag {mustBeTextScalar} = ""
%     options.resampling (1,1) double {mustBeMember(options.resampling,[0,1])} = 0
%     options.zeroPadding (1,1) double {mustBeInteger} = 0
%     options.conductivities (1,1) struct = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
%     options.voxSize (1,3) double = [1.0 1.0 1.0];
%     options.suppressMeshParameterWarning (1,1) logical = false;
%     options.suppressConductivityParameterWarning (1,1) logical = false;
%     options.customElectrodesTag (1,:) char = 'customLocations';
%     options.visualizeResult (1,1) logical = true;
%
% Description:
% Main function of ROAST.
% 
% Please refer to the README.md on the github repo for better formated
% documentations: https://github.com/andypotatohy/roast
% 
% If you use ROAST in your research, please cite these:
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., Realistic vOlumetric-Approach
% to Simulate Transcranial Electric Stimulation -- ROAST -- a fully automated
% open-source pipeline, Journal of Neural Engineering, Vol. 16, No. 5, 2019 (prefered reference)
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., ROAST: an open-source,
% fully-automated, Realistic vOlumetric-Approach-based Simulator for TES,
% Proceedings of the 40th Annual International Conference of the IEEE Engineering
% in Medicine and Biology Society, Honolulu, HI, July 2018
% 
% If you use New York head to run simulation, please also cite the following:
% Huang, Y., Parra, L.C., Haufe, S.,2016. The New York Head - A precise
% standardized volume conductor model for EEG source localization and tES
% targeting, NeuroImage,140, 150-162
% 
% If you also use the targeting feature (`roast_target`), please cite these:
% 
% Dmochowski, J.P., Datta, A., Bikson, M., Su, Y., Parra, L.C., Optimized 
% multi-electrode stimulation increases focality and intensity at target,
% Journal of Neural Engineering 8 (4), 046011, 2011
% 
% Dmochowski, J.P., Datta, A., Huang, Y., Richardson, J.D., Bikson, M.,
% Fridriksson, J., Parra, L.C., Targeted transcranial direct current stimulation 
% for rehabilitation after stroke, NeuroImage, 75, 12-19, 2013
% 
% Huang, Y., Thomas, C., Datta, A., Parra, L.C., Optimized tDCS for Targeting
% Multiple Brain Regions: An Integrated Implementation. Proceedings of the 40th
% Annual International Conference of the IEEE Engineering in Medicine and Biology
% Society, Honolulu, HI, July 2018, 3545-3548
% 
% ROAST was supported by NIH through grants R01MH111896, R01MH111439, 
% R01NS095123, R44NS092144, R41NS076123, and by Soterix Medical Inc.
% 
% General Public License version 3 or later. See LICENSE.md for details.
% 
% This software uses free packages from the Internet, except Matlab, which
% is a proprietary software by the MathWorks. You need a valid Matlab license
% to run this software.
% 
% ROAST is considered as an "aggregate" rather than "derived work", based on
% the definitions in GPL FAQ. The ROAST license only applies to the scripts,
% documentation and the individual MRI data under example/ folder in this 
% package and excludes those programs stored in the lib/ directory. The software 
% under lib/ follow their respective licenses. This software is only intended
% for non-commercial use.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% September 2019

arguments
    subj {mustBeTextScalar} = 'example/MNI152_T1_1mm.nii';
    recipe (1,:) cell = {'Fp1',1,'P4',-1};
    options.capType {mustBeMember(options.capType, {'1020','1010','1005','biosemi','egi'})} = '1010';
    options.elecType = [];
    options.elecSize = [];
    options.elecOri = 'lr'; % Can be 'lr', 'ap', 'si', or (1,3) vector giving precise vector of long axis
    options.T2 = [];
    options.meshOptions (1,1) struct = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);
    options.simulationTag {mustBeTextScalar} = ""
    options.resampling (1,1) double {mustBeMember(options.resampling,[0,1])} = 0
    options.zeroPadding (1,1) double {mustBeInteger} = 0
    options.conductivities (1,1) struct = struct('white', 0.126, 'gray', 0.276, 'csf', 1.65, 'bone', 0.01, 'skin', 0.465, 'air', 2.5e-14, 'gel', 0.3, 'electrode', 5.9e7);
    options.voxSize (1,3) double = [1.0 1.0 1.0];
    options.suppressMeshParameterWarning (1,1) logical = false;
    options.suppressConductivityParameterWarning (1,1) logical = false;
    options.customElectrodesTag (1,:) char = 'customLocations';
    options.visualizeResult (1,1) logical = true;
end

fprintf('\n\n');
disp('=============================================================')
disp('ROAST is an aggregated work by Yu (Andy) Huang licensed under')
disp('General Public License version 3 or later. It''s supported by')
disp('both NIH grants and Soterix Medical Inc.')
disp('=============================================================')

addpath(genpath([fileparts(which(mfilename)) filesep 'lib/']));

fprintf('\n\n');
disp('======================================================')
disp('CHECKING INPUTS...')
disp('======================================================')
fprintf('\n');

% check subject name
if nargin<1 || isempty(subj)
    subj = 'example/MNI152_T1_1mm.nii';
end

if strcmpi(subj,'nyhead')
    subj = 'example/nyhead.nii';
elseif strcmpi(subj, 'forrest')
    subj = 'example/forrest.nii';
elseif strcmpi(subj, 'forrestLarge')
    subj = 'example/forrestLarge.nii';
elseif ~exist(subj,'file')
    error('The subject MRI you provided ("%s") does not exist.', subj);
end

% Pull into own variables for further custom-error-checking:   
capType = options.capType;
elecType = options.elecType;
elecSize = options.elecSize;
elecOri = options.elecOri;
meshOptions = options.meshOptions;
if strlength(options.simulationTag) == 0
    uniqueTag = [];
else
    uniqueTag = options.simulationTag;
end
resampling = options.resampling;
zeroPadding = options.zeroPadding;
conductivities = options.conductivities;
voxSize = options.voxSize;

[elecPara, elecName, injectCurrent] = ...
    parse_electrode_parameters(recipe, capType, elecType, elecSize, elecOri);

if ~isempty(options.T2)
    if isstring(options.T2) || ischar(options.T2)
        if exist(options.T2,'file')==0
            error("The T2 MRI filename you provided (%s) does not exist.", options.T2);
        else
            T2 = load_untouch_nii(options.T2);
            if T2.hdr.hist.qoffset_x == 0 && T2.hdr.hist.srow_x(4)==0
                error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM Display function to fix the header.');
            end
        end
    else
        T2 = options.T2;
    end
else
    T2 = options.T2;
end
   

meshOptions = error_check_mesh_options(meshOptions, ...
    'SuppressWarning', options.suppressMeshParameterWarning);
conductivities = error_check_conductivities(conductivities, numel(elecName), ...
    'SuppressWarning', options.suppressConductivityParameterWarning);

% preprocess MRI data
if ~strcmpi(subj,'example/nyhead.nii') && ~strcmpi(subj, 'example/forrest.nii') && ~strcmpi(subj, 'example/forrestLarge.nii') % only when it's not NY head
    
    t1Data = load_untouch_nii(subj);
    if t1Data.hdr.hist.qoffset_x == 0 && t1Data.hdr.hist.srow_x(4)==0
        error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM Display function to fix the header.');
    end
    % check if bad MRI header

    if any(t1Data.hdr.dime.pixdim(2:4)<0.8) && ~resampling
        warning('The MRI has higher resolution (<0.8mm) in at least one direction. This will make the modeling process more computationally expensive and thus slower. If you wish to run faster using just 1-mm model, you can ask ROAST to re-sample the MRI into 1 mm first, by turning on the ''resampling'' option.');
    end
    % check if high-resolution MRI (< 0.8 mm in any direction)
    
    if length(unique(t1Data.hdr.dime.pixdim(2:4)))>1 && ~resampling
        warning('The MRI has anisotropic resolution. It is highly recommended that you turn on the ''resampling'' option, as the electrode size will not be exact if the model is built from an MRI with anisotropic resolution.');
    end
    % check if anisotropic resolution MRI
    
    [subjRas,isNonRAS] = convertToRAS(subj);
    % check if in non-RAS orientation, and if yes, put it into RAS
    
%     [subjRasRS,doResamp] = resampToOneMM(subjRas,doResamp);    
    [subjRasRS,resampling] = resampMesh(subjRas,resampling,voxSize);
    
    if zeroPadding>0
        subjRasRSPD = zeroPadding(subjRasRS,zeroPadding);
    else
        subjRasRSPD = subjRasRS;
    end
    
    if ~isempty(T2)
        T2 = realignT2(T2,subjRasRSPD);
    end
    % check if T2 is aligned with T1
    
else
    if strcmpi(subj, 'example/nyhead.nii')
        if ~exist('example/nyhead_T1orT2_masks.nii','file')
            unzip('example/nyhead_T1orT2_masks.nii.zip','example')
        end
        
        isNonRAS = 0; % New York head is in RAS
        
        if resampling
            error('The beauty of New York head is its 0.5 mm resolution. It''s a bad practice to resample it into 1 mm. Use another head ''example/MNI152_T1_1mm.nii'' for 1 mm model.');
        end
        
        if zeroPadding>0
            zeroPadding('example/nyhead_T1orT2_masks.nii',zeroPadding);
            subjRasRSPD = ['example/nyhead_padded' num2str(zeroPadding) '.nii'];
            if ~exist(['example/nyhead_padded' num2str(zeroPadding) '_T1orT2_seg8.mat'],'file')
                load('example/nyhead_T1orT2_seg8.mat','image','tpm','Affine');
                origin = inv(image.mat)*[0;0;0;1];
                origin = origin(1:3) + zeroPadding;
                image.mat(1:3,4) = [-dot(origin,image.mat(1,1:3));-dot(origin,image.mat(2,1:3));-dot(origin,image.mat(3,1:3))];
                save(['example/nyhead_padded' num2str(zeroPadding) '_T1orT2_seg8.mat'],'image','tpm','Affine');
            end
        else
            subjRasRSPD = subj;
        end
        
        if ~isempty(T2)
           warning('New York head selected. Any specified T2 image will be ignored.');
           T2 = [];
        end
    else
        isNonRAS = 0;
        if zeroPadding>0
            zeroPadding('example/forrest.nii',zeroPadding);
            subjRasRSPD = ['example/forrest_padded' num2str(zeroPadding) '.nii'];
            if ~exist(['example/forrest_padded' num2str(zeroPadding) '_T1orT2_seg8.mat'],'file')
                load('example/forrest_T1orT2_seg8.mat','image','tpm','Affine');
                origin = inv(image.mat)*[0;0;0;1];
                origin = origin(1:3) + zeroPadding;
                image.mat(1:3,4) = [-dot(origin,image.mat(1,1:3));-dot(origin,image.mat(2,1:3));-dot(origin,image.mat(3,1:3))];
                save(['example/forrest_padded' num2str(zeroPadding) '_T1orT2_seg8.mat'],'image','tpm','Affine');
            end
        else
            subjRasRSPD = subj;
        end
        T2 = [];
    end
        
end

% preprocess electrodes
[elecPara,indInUsrInput] = elecPreproc(subj,elecName,elecPara);

if any(~strcmpi(recipe,'leadfield'))
    
    elecName = elecName(indInUsrInput);
    injectCurrent = injectCurrent(indInUsrInput);
    
    configTxt = [];
    for i=1:length(elecName)
        configTxt = [configTxt elecName{i} ' (' num2str(injectCurrent(i)) ' mA), '];
    end
    configTxt = configTxt(1:end-2);
    
else
    
    elecNameOri = elecName; % back up for re-ordering solutions back to .loc file order;
                            % this is ugly, as .loc file has a different order of electrodes
                            % for historical reasons;
                            % HDE follows .loc file; ROAST follows capInfo.xls
    elecName = elecName(indInUsrInput);
    configTxt = 'leadFieldGeneration';
    
end

conductivities.gel = conductivities.gel(indInUsrInput);
conductivities.electrode = conductivities.electrode(indInUsrInput);

% sort elec options
if length(elecPara)==1
    if size(elecSize,1)>1, elecPara.elecSize = elecPara.elecSize(indInUsrInput,:); end
    if ~ischar(elecOri) && size(elecOri,1)>1
        elecPara.elecOri = elecPara.elecOri(indInUsrInput,:);
    end
elseif length(elecPara)==length(elecName)
    elecPara = elecPara(indInUsrInput);
else
    error('Something is wrong!');
end

opts = struct('configTxt',configTxt, ...
              'elecPara',elecPara, ...
              'T2',T2, ...
              'meshOpt',meshOptions, ...
              'conductivities',conductivities, ...
              'uniqueTag',uniqueTag, ...
              'resamp',resampling, ...
              'zeroPad',zeroPadding, ...
              'isNonRAS',isNonRAS, ...
              'customElectrodesTag',options.customElectrodesTag);

% log tracking
[dirname,baseFilename] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

Sopt = dir([dirname filesep baseFilename '_*_roastOptions.mat']);
if isempty(Sopt)
    opts = writeRoastLog(subj,opts,'roast');
else
    isNew = zeros(length(Sopt),1);
    for i=1:length(Sopt)
        load([dirname filesep Sopt(i).name],'opt');
        isNew(i) = isNewOptions(opts,opt,'roast');
    end
    if all(isNew)
        opts = writeRoastLog(subj,opts,'roast');
    else
        load([dirname filesep Sopt(find(~isNew)).name],'opt');
        if ~isempty(opts.uniqueTag) && ~strcmp(opts.uniqueTag,opt.uniqueTag)
            warning(['The simulation with the same options has been run before under tag ''' opt.uniqueTag '''. The new tag you specified ''' opts.uniqueTag ''' will be ignored.']);
        end
        opts.uniqueTag = opt.uniqueTag;
    end
end
% uniqueTag = opts.uniqueTag;

fprintf('\n');
disp('======================================================')
if strcmp(baseFilename,'nyhead')
    disp('ROAST New York head')
elseif strcmp(baseFilename, 'forrest')
    disp('ROAST Forrest NHP head');
elseif strcmp(baseFilename, 'forrestLarge')
    disp('ROAST Forrest LARGE ELECTRODE model');
else
    disp(['ROAST ' subj])
end
disp('USING RECIPE:')
disp(configTxt)
disp('...and simulation options saved in:')
disp([dirname filesep baseFilename '_roastLog,'])
disp(['under tag: ' uniqueTag])
disp('======================================================')
fprintf('\n\n');

% warn users lead field will take a long time to generate
if all(strcmpi(recipe,'leadfield'))
    [~,indRef] = ismember('Iz',elecName);
    indStimElec = setdiff(1:length(elecName),indRef);
    [isInRoastCore,indInRoastCore] = ismember(elecNameOri,elecName(indStimElec));
    isSolved = zeros(length(indStimElec),1);
    for i=1:length(indStimElec)
        if exist([dirname filesep baseFilename '_' uniqueTag '_e' num2str(indStimElec(i)) '.pos'],'file')
            isSolved(i) = 1;
        end
    end
    % only warn users the first time they run for this subject
    if all(~isSolved) && ~exist([dirname filesep baseFilename '_' uniqueTag '_roastResult.mat'],'file')
        warning('You specified the ''recipe'' as the ''lead field generation''. Nice choice! Note all customized options on electrodes are overwritten by the defaults. Refer to the readme file for more details. Also this will usually take a long time (>1 day) to generate the lead field for all the candidate electrodes.');
        doLFconfirm = input('Do you want to continue? ([Y]/N)','s');
        if strcmpi(doLFconfirm,'n'), disp('Aborted.'); return; end
    end
end

if ~(strcmp(baseFilename,'nyhead') || strcmp(baseFilename, 'forrest') || strcmp(baseFilename, 'forrestLarge'))
    
    [~,baseFilenameRasRSPD] = fileparts(subjRasRSPD);
    
    if (isempty(T2) && ~exist([dirname filesep 'c1' baseFilenameRasRSPD '_T1orT2.nii'],'file')) ||...
            (~isempty(T2) && ~exist([dirname filesep 'c1' baseFilenameRasRSPD '_T1andT2.nii'],'file'))
        disp('======================================================')
        disp('       STEP 1 (out of 6): SEGMENT THE MRI...          ')
        disp('======================================================')
        start_seg(subjRasRSPD,T2);
    else
        disp('======================================================')
        disp('          MRI ALREADY SEGMENTED, SKIP STEP 1          ')
        disp('======================================================')
    end
    
    if (isempty(T2) && ~exist([dirname filesep baseFilenameRasRSPD '_T1orT2_masks.nii'],'file')) ||...
            (~isempty(T2) && ~exist([dirname filesep baseFilenameRasRSPD '_T1andT2_masks.nii'],'file'))
        disp('======================================================')
        disp('     STEP 2 (out of 6): SEGMENTATION TOUCHUP...       ')
        disp('======================================================')
        segTouchup(subjRasRSPD,T2);
    else
        disp('======================================================')
        disp('    SEGMENTATION TOUCHUP ALREADY DONE, SKIP STEP 2    ')
        disp('======================================================')
    end
    
else
    
    disp('======================================================')
    disp(' PRE-SEGMENTED SELECTED, GOING TO STEP 3 DIRECTLY...  ')
    disp('======================================================')
%     warning('New York head is a 0.5 mm model so is more computationally expensive. Make sure you have a decent machine (>50GB memory) to run ROAST with New York head.')
    [~,baseFilenameRasRSPD] = fileparts(subjRasRSPD);
    
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '_mask_elec.nii'],'file')
    disp('======================================================')
    disp('      STEP 3 (out of 6): ELECTRODE PLACEMENT...       ')
    disp('======================================================')
    hdrInfo = electrodePlacement(subj,subjRasRSPD,T2,elecName,opts,uniqueTag);
else
    disp('======================================================')
    disp('         ELECTRODE ALREADY PLACED, SKIP STEP 3        ')
    disp('======================================================')
%     load([dirname filesep baseFilename '_' uniqueTag '_labelVol.mat'],'volume_elecLabel','volume_gelLabel');
    load([dirname filesep baseFilenameRasRSPD '_header.mat'],'hdrInfo');
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '.mat'],'file')
    disp('======================================================')
    disp('        STEP 4 (out of 6): MESH GENERATION...         ')
    disp('======================================================')
    [node,elem,face] = meshByIso2mesh(subj,subjRasRSPD,T2,meshOptions,hdrInfo,uniqueTag);
else
    disp('======================================================')
    disp('          MESH ALREADY GENERATED, SKIP STEP 4         ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '.mat'],'node','elem','face');
end

if any(~strcmpi(recipe,'leadfield'))
    
    if ~exist([dirname filesep baseFilename '_' uniqueTag '_v.pos'],'file')
        disp('======================================================')
        disp('       STEP 5 (out of 6): SOLVING THE MODEL...        ')
        disp('======================================================')
        prepareForGetDP(subj,node,elem,elecName,uniqueTag);
        indElecSolve = 1:length(elecName);
        solveByGetDP(subj,injectCurrent,conductivities,indElecSolve,uniqueTag,'');
    else
        disp('======================================================')
        disp('           MODEL ALREADY SOLVED, SKIP STEP 5          ')
        disp('======================================================')
        %     load([dirname filesep baseFilename '_' uniqueTag '_elecMeshLabels.mat'],'label_elec');
    end
    
    if ~exist([dirname filesep baseFilename '_' uniqueTag '_roastResult.mat'],'file')
        if options.visualizeResult
            disp('======================================================')
            disp('STEP 6 (final step): SAVING AND VISUALIZING RESULTS...')
            disp('======================================================')
            [vol_all,ef_mag,ef_all] = postGetDP(subj,subjRasRSPD,node,hdrInfo,uniqueTag);
            visualizeRes(subj,subjRasRSPD,T2,node,elem,face,injectCurrent,hdrInfo,uniqueTag,0,vol_all,ef_mag,ef_all);
        else
            disp('======================================================')
            disp('STEP 6 (final step):            FORCE-SKIPPED!        ')
            disp('======================================================')
        end
    else
        disp('======================================================')
        disp('  ALL STEPS DONE, LOADING RESULTS FOR VISUALIZATION   ')
        disp('======================================================')
        load([dirname filesep baseFilename '_' uniqueTag '_roastResult.mat'],'vol_all','ef_mag','ef_all');
        visualizeRes(subj,subjRasRSPD,T2,node,elem,face,injectCurrent,hdrInfo,uniqueTag,1,vol_all,ef_mag,ef_all);
    end
    
else
    
    if any(~isSolved) && ~exist([dirname filesep baseFilename '_' uniqueTag '_roastResult.mat'],'file')
        disp('======================================================')
        disp('    STEP 5 (out of 6): GENERATING THE LEAD FIELD...   ')
        disp('           NOTE THIS WILL TAKE SOME TIME...           ')
        disp('======================================================')
        prepareForGetDP(subj,node,elem,elecName,uniqueTag);
        injectCurrent = ones(length(elecName),1); % 1 mA at each candidate electrode
        injectCurrent(indRef) = -1;
        for i=1:length(indStimElec)
            if ~isSolved(i)
                fprintf('\n======================================================\n');
                disp(['SOLVING FOR ELECTRODE ' num2str(i) ' OUT OF ' num2str(length(indStimElec)) ' ...']);
                fprintf('======================================================\n\n');
                indElecSolve = [indStimElec(i) indRef];
                solveByGetDP(subj,injectCurrent,conductivities,indElecSolve,uniqueTag,num2str(indStimElec(i)));
            else
                disp(['ELECTRODE ' num2str(i) ' HAS BEEN SOLVED, SKIPPING...']);
            end
        end
    else
        disp('======================================================')
        disp('       LEAD FIELD ALREADY GENERATED, SKIP STEP 5      ')
        disp('======================================================')
        %     load([dirname filesep baseFilename '_' uniqueTag '_elecMeshLabels.mat'],'label_elec');
    end
    
    if ~exist([dirname filesep baseFilename '_' uniqueTag '_roastResult.mat'],'file')
        disp('========================================================')
        disp('STEP 6 (final step): ASSEMBLING AND SAVING LEAD FIELD...')
        disp('========================================================')
        postGetDP(subj,[],node,hdrInfo,uniqueTag,indStimElec,indInRoastCore(isInRoastCore));
    else
        disp('======================================================')
        disp('         ALL STEPS DONE, READY TO DO TARGETING        ')
        disp(['         FOR SUBJECT ' subj])
        disp(['         USING TAG ' uniqueTag])
        disp('======================================================')
    end
    
end

disp('==================ALL DONE ROAST=======================');