% Run the LEADS MRI-based PET processing pipeline with user input on the
% tasks to be completed. Full pipeline includes:
%
% Setup
% - Newdata files unzipped
% - DICOMs converted to NIfTI
% - MRI and PET scans moved from newdata to raw directories
% - CSV files saved indicating which scans are scheduled for processing
% - Processed PET and MRI directories created and linked appropriately
%
% T1 MRI processing:
% - MRIs processed through FreeSurfer 7.1 recon-all + brainstem seg
% - MRIs and FreeSurfer-based parcellation files recentered to the
%   nu.nii center-of-mass then coregistered to SPM12 OldNorm T1 template
% - Follow-up MRIs coregistered to baseline MRI
% - MRIs segmented and warped to MNI space
% - ROI masks saved out for reference regions and target ROIs
% - MRIs affine transformed to transformation to the SPM12 OldNorm T1
%   template (for visualization purposes only)
%
% PET processing:
% - PET scans for all 3 tracers (FBB, FDG, FTP) are processed using the
%   same steps, only different reference regions. This part of the
%   processing must happen after MRI processing is complete.
% - PET scans recentered to axis midpoints, then coregistered to the
%   closest MRI (you can always check which MRI was used by looking at
%   the 'mri' symlink in each processed PET directory)
% - PET scans intensity normalized to tracer-specific reference regions
%   based on FreeSurfer regions to make voxelwise SUVR images. In some
%   cases, multiple SUVRs may be saved for a single PET scan.
%     - FBB: whole cerebellum, composite white matter (ADNI method)
%     - FDG: pons
%     - FTP: inferior cerebellar gray matter, eroded subcortical white
%            matter (Berkeley method)
% - PET SUVRs warped to MNI space using forward deformation files from
%   MRI segmentation
% - PET SUVRs affine transformed to MNI space using transform estimated
%   for the nu.nii (for visualization purposes only)
%
% ROI extraction:
% - PET SUVR means calculated for all FreeSurfer regions in the
%   aparc+aseg.nii
%
%                                  Processing code by D. Schonhaut, 2024
% ----------------------------------------------------------------------

% Start the program timer
tic;

% Say hi
say_hi();

% Reset the MATLAB environment
restoredefaultpath;
addpath("/home/mac/dschonhaut/code/matlab/spm12");
addpath("/home/mac/dschonhaut/code/matlab/spm12/toolbox");
addpath("/home/mac/dschonhaut/code/matlab/spm12/toolbox/suit");
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

% Define defaults
proj_dir = "/mnt/coredata/processing/leads";
code_dir = fileparts(fileparts(mfilename('fullpath')));
scans_to_process_dir = fullfile(proj_dir, 'metadata', 'scans_to_process');
skip_newdata = false;
wipe_newdata = true;
setup_processed = true;
dry_run = false;
schedule_all_mris = false;
overwrite_raw = false;
schedule_overwrite = false;
overwrite_processed = false;
overwrite = false;
mri_dirs = {};
segment_brainstem = true;
process_freesurfer = true;
process_post_freesurfer = true;
pet_dirs = {};

% Get path to the Python interpreter and Python scripts that this
% program can call
python = '/home/mac/dschonhaut/mambaforge/envs/nipy311/bin/python';
qc_evals_script = fullfile(code_dir, 'qc', 'qc_evals.py');

% Figure out what the user wants to do
prompt_user = sprintf([
    '\n~ Welcome to the LEADS processing pipeline ~\n\n', ...
    'What would you like to do?\n', ...
    '  [1] Add scans from newdata to raw, and\n', ...
    '      schedule raw scans to be processed\n', ...
    '  [2] View scheduled MRI and PET scans\n', ...
    '      (but don''t process them yet)\n', ...
    '  [3] Process scheduled MRIs\n', ...
    '  [4] Process scheduled PET scans\n', ...
    '  [5] Merge QC evaluation files\n', ...
    '  [6] View processed scans that need to be QC''d\n', ...
    '  [7] Merge ROI extraction files (NOT IMPLEMENTED) \n', ...
    '  [8] Prepare quarterly report files (NOT IMPLEMENTED) \n', ...
    '  [9] Exit\n\n'
]);
action = input(prompt_user);

% Run the selected action
change_defaults_msg = '\nDo you need to change any defaults?';
overwrite_msg = 'Overwrite existing files?';
switch action
    case 1
        % Schedule scans to be processed
        fprintf('\nDefault parameters\n------------------\n');
        fprintf('Project directory                                                  = %s\n', proj_dir);
        fprintf('Skip the ''newdata'' -> ''raw'' submodule?                             = %d\n', skip_newdata);
        fprintf('Wipe ''newdata'' after moving scans to ''raw''                         = %d\n', wipe_newdata);
        fprintf('Do a dry run (summarize but do not schedule scans yet)             = %d\n', dry_run);
        fprintf('Schedule all unprocessed MRIs in ''raw'' (including orphans)         = %d\n', schedule_all_mris);
        fprintf('Setup ''processed'' directories for scheduled scans                  = %d\n', setup_processed);
        fprintf('Overwrite existing directories when moving ''newdata'' -> ''raw''      = %d\n', overwrite_raw);
        fprintf('Schedule already processed scans to be reprocessed                 = %d\n', schedule_overwrite);
        fprintf('Remove processed directories for scans scheduled to be reprocessed = %d\n', overwrite_processed);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            require_response = false;
            confirm_response = true;
            proj_dir = prompt_text('Enter path to project directory', proj_dir, confirm_response);
            skip_newdata = prompt_bool('Skip ''newdata'' -> ''raw''?', false);
            if ~skip_newdata
                wipe_newdata = prompt_bool('Wipe ''newdata'' after moving scans?', true);
            end
            dry_run = prompt_bool('Do a dry run?', dry_run);
            schedule_all_mris = prompt_bool('Schedule all unprocessed MRIs?', schedule_all_mris);
            if ~dry_run
                setup_processed = prompt_bool('Setup scheduled scans in ''processed''?', true);
                if ~skip_newdata
                    overwrite_raw = prompt_bool('Overwrite scans in ''raw''?', overwrite_raw, require_response, confirm_response);
                end
                schedule_overwrite = prompt_bool('Schedule scans to be reprocessed?', schedule_overwrite, require_response, confirm_response);
                if setup_processed
                    overwrite_processed = prompt_bool('Overwrite scans in ''processed''?', overwrite_processed, require_response, confirm_response);
                end
            end
        end
        proj_dir = abspath(proj_dir);
        mustBeFolder(proj_dir);
        if schedule_overwrite && setup_processed && overwrite_processed
            proceed_knowingly = prompt_bool('!! PLEASE BE CAREFUL. You are about to remove and reset already processed scans in ''processed''. Are you sure you want to do this?', false, true, true);
            if ~proceed_knowingly
                fprintf('OK--let''s get out of this scary situation...\n');
                return;
            end
        end
        fprintf('\nCalling setup_leads_processing.m...\n');
        setup_leads_processing( ...
            proj_dir, ...
            skip_newdata, ...
            wipe_newdata, ...
            setup_processed, ...
            dry_run, ...
            schedule_all_mris, ...
            overwrite_raw, ...
            schedule_overwrite, ...
            overwrite_processed ...
        );
    case 2
        % View scheduled MRI and PET scans
        fprintf('\nCalling queue_mris_to_process.m...\n');
        queue_mris_to_process(scans_to_process_dir);
        fprintf('\nCalling queue_pets_to_process.m...\n');
        queue_pets_to_process(scans_to_process_dir);
    case 3
        % Process scheduled MRIs
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('MRIs to process                                   = By default, this program will attempt to process all    \n');
        fprintf('                                                    scheduled MRIs, but it is possible to override this     \n');
        fprintf('                                                    behavior and provide an exact list of scans you want to \n');
        fprintf('                                                    process. Note that if you intend to overwrite already   \n');
        fprintf('                                                    processed files, you still need to specify overwrite = 1\n');
        fprintf('                                                    or nothing will happen                                  \n');
        fprintf('Overwrite existing files                          = %d\n', overwrite);
        fprintf('Directory with CSV files listing scans to process = %s\n', scans_to_process_dir);
        fprintf('Run FreeSurfer processing                         = %d\n', process_freesurfer);
        fprintf('Segment brainstem                                 = %d\n', segment_brainstem);
        fprintf('Run post-FreeSurfer, SPM-based processing         = %d\n', process_post_freesurfer);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            set_mri_dirs = prompt_bool('Specify a list of MRIs to process?', false);
            if set_mri_dirs
                response = prompt_text( ...
                    'Enter path to the first MRI you want to process', ...
                    '', ...
                    false ...
                );
                if ~isfolder(response)
                    fprintf('\n!! WARNING: %s is not a folder and will not be added to the list\n\n', response);
                else
                    mri_dirs = [mri_dirs, response]
                end

                while 1
                    msg = sprintf('Enter path to the next MRI you want to process, or type ''q'' to move on\n');
                    response = prompt_text(msg, '', false);
                    if strcmp(lower(response(1)), 'q')
                        break;
                    elseif ~isfolder(response)
                        fprintf('\n!! WARNING: %s is not a folder and will not be added to the list\n\n', response);
                    else
                        mri_dirs = [mri_dirs, response]
                    end
                end
                mri_dirs = unique(abspath(cellvec(mri_dirs)));

                % Confirm the submitted directories
                n_scans = length(mri_dirs);
                fprintf('\nYou have submitted the following %d MRI scan directories to process:\n', n_scans);
                for i = 1:n_scans
                    fprintf('  %s\n', mri_dirs{i});
                end
                confirm_response = prompt_bool('Is this correct?', false, true);
                if ~confirm_response
                    fprintf('OK--let''s exit the program and try again...\n');
                    return;
                end
            end
            overwrite = prompt_bool(overwrite_msg, false, false, true);
            if isempty(mri_dirs)
                scans_to_process_dir = prompt_text('Enter path to ''scans_to_process'' directory', scans_to_process_dir, false);
            end
            process_freesurfer = prompt_bool('Process MRIs through FreeSurfer?', true);
            if process_freesurfer
                segment_brainstem = prompt_bool('Segment brainstem?', true);
            end
            process_post_freesurfer = prompt_bool('Run post-FreeSurfer MRI processing?', true);
        end
        fprintf('\nCalling process_mris.m...\n');
        process_mris( ...
            mri_dirs, ...
            scans_to_process_dir, ...
            overwrite, ...
            segment_brainstem, ...
            process_freesurfer, ...
            process_post_freesurfer ...
        );
    case 4
        % Process scheduled PET scans
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('PET scans to process                              = By default, this program will attempt to process all     \n');
        fprintf('                                                    scheduled PETs scans, but it is possible to override this\n');
        fprintf('                                                    behavior and provide an exact list of scans you want to  \n');
        fprintf('                                                    process. Note that if you intend to overwrite already    \n');
        fprintf('                                                    processed files, you still need to specify overwrite = 1 \n');
        fprintf('                                                    or nothing will happen                                   \n');
        fprintf('Overwrite existing files                          = %d\n', overwrite);
        fprintf('Directory with CSV files listing scans to process = %s\n', scans_to_process_dir);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            set_pet_dirs = prompt_bool('Specify a list of PET scans to process?', false);
            if set_pet_dirs
                response = prompt_text( ...
                    'Enter path to the first PET scan you want to process', ...
                    '', ...
                    false ...
                );
                if ~isfolder(response)
                    fprintf('\n!! WARNING: %s is not a folder and will not be added to the list\n\n', response);
                else
                    pet_dirs = [pet_dirs, response]
                end

                while 1
                    msg = sprintf('Enter path to the next PET scan you want to process, or type ''q'' to move on\n');
                    response = prompt_text(msg, '', false);
                    if strcmp(lower(response(1)), 'q')
                        break;
                    elseif ~isfolder(response)
                        fprintf('\n!! WARNING: %s is not a folder and will not be added to the list\n\n', response);
                    else
                        pet_dirs = [pet_dirs, response]
                    end
                end
                pet_dirs = unique(abspath(cellvec(pet_dirs)));

                % Confirm the submitted directories
                n_scans = length(pet_dirs);
                fprintf('\nYou have submitted the following %d PET scan directories to process:\n', n_scans);
                for i = 1:n_scans
                    fprintf('  %s\n', pet_dirs{i});
                end
                confirm_response = prompt_bool('Is this correct?', false, true);
                if ~confirm_response
                    fprintf('OK--let''s exit the program and try again...\n');
                    return;
                end
            end
            overwrite = prompt_bool(overwrite_msg, false, false, true);
            if isempty(pet_dirs)
                scans_to_process_dir = prompt_text('Enter path to ''scans_to_process'' directory', scans_to_process_dir, false);
            end
        end
        fprintf('\nCalling process_pets.m...\n');
        process_pets( ...
            pet_dirs, ...
            scans_to_process_dir, ...
            overwrite ...
        );
    case 5
        % Merge QC evaluation files
        cmd = sprintf('%s %s %s', python, qc_evals_script, 'merge');
        fprintf('Merging QC evaluation files...\n');
        fprintf('$ %s\n', cmd);
        system(cmd);
    case 6
        % View processed scans that need to be QC'd
        fprintf([ ...
            '\n** NOTE: Merge QC evaluation files (Option 5) first to view an\n', ...
            '         updated list of processed scans that need to be QC''d\n\n' ...
        ]);
        cmd = sprintf('%s %s %s', python, qc_evals_script, 'incomplete');
        fprintf('$ %s\n', cmd);
        system(cmd);
    case 7
        fprintf('Option not yet implemented. Exiting program.\n');
    case 8
        fprintf('Option not yet implemented. Exiting program.\n');
    case 9
        fprintf('Bye for now\n');
    otherwise
        fprintf('Input is invalid. Exiting program.\n');
end

% End the program
say_bye();
tock();
fprintf('\n');
clear;
