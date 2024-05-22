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
% - ROI masks saved out for reference regions and target ROIs
% - MRIs segmented and warped to MNI space
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
addpath(genpath("/mnt/coredata/processing/leads/code"));

% Define defaults
proj_dir = "/mnt/coredata/processing/leads";
scans_to_process_dir = fullfile(proj_dir, 'metadata', 'scans_to_process');
overwrite = false;
cleanup_newdata = true;
process_all_mris = false;
segment_brainstem = true;
process_freesurfer = true;
process_post_freesurfer = true;

% Figure out what the user wants to do
prompt_user = sprintf([
    '\n~ Welcome to the LEADS processing pipeline ~\n\n', ...
    'What do you want to do?\n', ...
    '  [1] Setup scans for processing\n', ...
    '  [2] View scans that are scheduled for processing (but don''t do anything)\n', ...
    '  [3] Process MRIs\n', ...
    '  [4] Process PET scans\n', ...
    '  [5] Nothing, I was just curious what happens when I run this script\n', ...
    '>> '
]);
action = input(prompt_user);

% Run the selected action
change_defaults_msg = 'Do you need to change any defaults?';
overwrite_msg = 'Overwrite existing files?';
switch action
    case 1
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('  Project directory                           : %s\n', proj_dir);
        fprintf('  Overwrite existing files                    : %d\n', overwrite);
        fprintf('  Cleanup newdata after moving scans          : %d\n', cleanup_newdata);
        fprintf('  Process all MRIs in raw (including orphans) : %d\n', process_all_mris);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            proj_dir = prompt_text('Enter path to project directory', proj_dir, true);
            overwrite = prompt_bool(overwrite_msg, false, false, true);
            cleanup_newdata = prompt_bool('Wipe ''newdata'' after moving scans to ''raw''?', true);
            process_all_mris = prompt_bool('Process all MRIs in ''raw'' (even orphans)?', false);
        end
        proj_dir = abspath(proj_dir);
        mustBeFolder(proj_dir);
        setup_leads_processing(proj_dir, overwrite, cleanup_newdata, process_all_mris);
    case 2
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('  Directory with CSV files listing scans to process : %s\n', scans_to_process_dir);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            scans_to_process_dir = prompt_text('Enter path to ''scans_to_process'' directory', scans_to_process_dir, false);
        end
        queue_mris_to_process(scans_to_process_dir);
        queue_pets_to_process(scans_to_process_dir);
    case 3
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('  Overwrite existing files                          : %d\n', overwrite);
        fprintf('  Directory with CSV files listing scans to process : %s\n', scans_to_process_dir);
        fprintf('  Segment brainstem                                 : %d\n', segment_brainstem);
        fprintf('  Run FreeSurfer processing                         : %d\n', process_freesurfer);
        fprintf('  Run post-FreeSurfer processing                    : %d\n', process_post_freesurfer);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            overwrite = prompt_bool(overwrite_msg, false, false, true);
            scans_to_process_dir = prompt_text('Enter path to ''scans_to_process'' directory', scans_to_process_dir, false);
            segment_brainstem = prompt_bool('Segment brainstem?', true);
            process_freesurfer = prompt_bool('Process MRIs through FreeSurfer?', true);
            process_post_freesurfer = prompt_bool('Run post-FreeSurfer MRI processing?', true);
        end
        process_mris(overwrite, scans_to_process_dir, segment_brainstem, process_freesurfer, process_post_freesurfer);
    case 4
        fprintf('\nDefaults parameters\n-------------------\n');
        fprintf('  Overwrite existing files                          : %d\n', overwrite);
        fprintf('  Directory with CSV files listing scans to process : %s\n', scans_to_process_dir);
        change_defaults = prompt_bool(change_defaults_msg, false);
        if change_defaults
            overwrite = prompt_bool(overwrite_msg, false, false, true);
            scans_to_process_dir = prompt_text('Enter path to ''scans_to_process'' directory', scans_to_process_dir, false);
        end
        process_pets(overwrite, scans_to_process_dir);
    case 5
        fprintf('Ok\n');
    otherwise
        fprintf('INVALID INPUT. Try again :(\n');
end

% End the program
say_bye();
tock();
clear;
