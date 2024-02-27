%%%%% Master script to process LEADS PET data %%%%%
%%%%% Original code by Leo Iaccarino in 2019 %%%%%
%%%%% Updated by Daniel Schonhaut in 2024 %%%%%

function removeSubdirs(paths, prefix, verbose)
% Remove subdirectories with a specified prefix in given paths.
%
% Parameters
% ----------
%   paths : cell array of strings
%     Each string is a path in which directories with the specified
%     prefix will be searched for and deleted. Required.
%   prefix : string, optional
%     String specifying the prefix of the directories to delete.
%     Default is '' (all subdirectories are deleted).
%   verbose : logical, optional
%     If true, prints the path of each directory that is deleted.
%     Default is true.
%

% Set default values for optional parameters.
switch nargin
    case 1
        prefix = '';
        verbose = true;
    case 2
        verbose = true;
end

% Ensure paths is a cell array.
if ischar(paths) || isstring(paths)
    paths = {paths};
end

% Remove subdirectories with the specified prefix.
for ii = 1:length(paths)
    parentDir = paths{ii};
    if ~isfolder(parentDir)
        if verbose
            fprintf('Deleting directories; %s is not an existing directory, moving on...\n', parentDir);
        end
        continue;
    end
    subDirs = dir(fullfile(parentDir, [prefix '*']));
    for jj = 1:length(subDirs)
        if subDirs(jj).isdir
            subDir = fullfile(parentDir, subDirs(jj).name);
            % rmdir(subDir, 's');
            if verbose
                fprintf('Deleted %s\n', subDir);
            end
        end
    end
end
end


% Define directory paths
dirs = containers.Map;
dirs('proj') = '/mnt/coredata/processing/leads';
dirs('code') = fullfile(dirs('proj'), 'code');
dirs('data') = fullfile(dirs('proj'), 'data');
dirs('extraction') = fullfile(dirs('data'), 'extraction');
dirs('freesurfer') = fullfile(dirs('data'), 'freesurfer');
dirs('links') = fullfile(dirs('data'), 'links');
dirs('newdata') = fullfile(dirs('data'), 'newdata');
dirs('processed') = fullfile(dirs('data'), 'processed');

% Add paths to other scripts we need to access
addpath(genpath(dirs('code')));
addpath(genpath('/mnt/coredata/Projects/Resources/scripts/leotools'));

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
    ], dirs('proj')
);
user_action = input(msg);
if user_action == 1
    run LEADS_convert_dicoms.m
    run LEADS_MRI_Processing.m
    removeSubdirs(fullfile(dirs('newdata'), {'mri'}));
elseif user_action == 2
    run LEADS_convert_dicoms.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(dirs('newdata'), {'fbb', 'fdg', 'ftp'}));
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
    removeSubdirs(fullfile(dirs('newdata'), {'mri', 'fbb', 'fdg', 'ftp'}));
elseif user_action == 7
    run LEADS_convert_dicoms.m
    run LEADS_MRI_Processing.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(dirs('newdata'), {'mri', 'fbb', 'fdg', 'ftp'}));
    run LEADS_PET_Quantification.m
    run LEADS_PETMRI_MNI.m
    run LEADS_Wmapping.m
elseif user_action == 8
    run LEADS_convert_dicoms.m
    run LEADS_PET_Processing.m
    removeSubdirs(fullfile(dirs('newdata'), {'fbb', 'fdg', 'ftp'}));
    run LEADS_PET_Quantification.m
    run LEADS_PETMRI_MNI.m
    run LEADS_Wmapping.m
end

% End the program
msg = [
    '\nAll done',...
    '\n',...
    '\n   |\\      _,,,---,,_      ',...
    '\n   /,`.-''`''    -.  ;-;;,_  ',...
    '\n  |,4-  ) )-,_. ,\\ (  `''-'' ',...
    '\n ''---''''(_/--''  `-''\\_)      ',...
    '\n\n'
    ]
fprintf(msg);
clear;
