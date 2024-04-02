%%%%% Master script to process LEADS PET data %%%%%
%%%%% Original code by Leo Iaccarino in 2019 %%%%%
%%%%% Updated by Daniel Schonhaut in 2024 %%%%%

% Add paths to other scripts we need to access
restoredefaultpath;
run '/home/mac/dschonhaut/matlab/startup.m';

% Define directory paths
PATHS = containers.Map;
PATHS('proj') = '/mnt/coredata/processing/leads';
PATHS('code') = fullfile(PATHS('proj'), 'code');
PATHS('data') = fullfile(PATHS('proj'), 'data');
PATHS('extraction') = fullfile(PATHS('data'), 'extraction');
PATHS('freesurfer') = fullfile(PATHS('data'), 'freesurfer');
PATHS('links') = fullfile(PATHS('data'), 'links');
PATHS('metadata') = fullfile(PATHS('proj'), 'metadata');
PATHS('raw') = fullfile(PATHS('data'), 'raw');
PATHS('processed') = fullfile(PATHS('data'), 'processed');
addpath(genpath(PATHS('code')));

% Ask the user what processing to perform
msg = sprintf(
[
    '\nWelcome to the LEADS processing pipeline!',...
    '\n\nThe top-level directory is set to: %s',...
    '\n\nWhat do you want to do?',...
    '\n  [1] Process MRIs through FreeSurfer',...
    '\n  [2] Process PET through MRI-based pipeline',...
    '\n      (assumes MRIs have already been processed)',...
    '\n  [3] Extract PET SUVR means from FreeSurfer ROIs',...
    '\n  [4] Warp MRI and PET to MNI space',...
    '\n  [5] Create MRI and PET W-score maps',...
    '\n  [6] Process MRIs and then PET (actions [1] and [2])',...
    '\n  [7] Run the full pipeline (actions [1]-[6])',...
    '\n  [8] Run the full pipeline except MRI processing (actions [2]-[6])',...
    '\n  --> '
    ], PATHS('proj')
);
user_action = input(msg);
if user_action == 1
    run LEADS_convert_dicoms.m
    run LEADS_MRI_Processing.m
    removeSubdirs(fullfile(PATHS('newdata'), {'mri'}));
elseif user_action == 2
    run LEADS_convert_dicoms.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(PATHS('newdata'), {'fbb', 'fdg', 'ftp'}));
elseif user_action == 3
    run LEADS_PET_Quantification.m
elseif user_action == 4
    run LEADS_PETMRI_MNI.m
elseif user_action == 5
    run LEADS_Wmapping.m
elseif user_action == 6
    run LEADS_convert_dicoms.m
    run LEADS_MRI_Processing.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(PATHS('newdata'), {'mri', 'fbb', 'fdg', 'ftp'}));
elseif user_action == 7
    run LEADS_convert_dicoms.m
    run LEADS_MRI_Processing.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(PATHS('newdata'), {'mri', 'fbb', 'fdg', 'ftp'}));
    run LEADS_PET_Quantification.m
    run LEADS_PETMRI_MNI.m
    run LEADS_Wmapping.m
elseif user_action == 8
    run LEADS_convert_dicoms.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(PATHS('newdata'), {'fbb', 'fdg', 'ftp'}));
    run LEADS_PET_Quantification.m
    run LEADS_PETMRI_MNI.m
    run LEADS_Wmapping.m
end

% End the program
sayBye();
clear;
