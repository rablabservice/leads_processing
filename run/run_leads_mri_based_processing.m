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

% Reset the MATLAB environment
restoredefaultpath;
addpath("/home/mac/dschonhaut/code/matlab/spm12");
addpath("/home/mac/dschonhaut/code/matlab/spm12/toolbox");
addpath("/home/mac/dschonhaut/code/matlab/spm12/toolbox/suit");
addpath(genpath("/mnt/coredata/processing/leads/code"));

% Define globals
data_dir = "/mnt/coredata/processing/leads/data";
overwrite = false;
process_unused_mris = false;
cleanup = true;
log_dir = "/mnt/coredata/processing/leads/metadata/scans_to_process";
segment_brainstem = true;

% Say hi
tic;
say_hi();

% Figure out what the user wants to do
prompt_user = sprintf([
    '\n~ Welcome to the LEADS processing pipeline ~\n\n', ...
    'What do you want to do?\n', ...
    '  [1] Setup scans for processing\n', ...
    '  [2] Process MRIs\n', ...
    '  [3] Process PET scans\n', ...
    '  [4] Process MRIs and PET  ([2]-[3])\n', ...
    '  [5] Run the full pipeline ([1]-[3])\n', ...
    '  [6] Nothing, I was just curious what happens when I run this script\n', ...
    '>> '
]);
action = input(prompt_user);

% Run the selected action
switch action
    case 1
        setup_leads_processing(data_dir, overwrite, process_unused_mris, cleanup);
    case 2
        process_mris(overwrite, log_dir, segment_brainstem);
    case 3
        process_pets(overwrite, log_dir);
    case 4
        process_mris(overwrite, log_dir, segment_brainstem);
        process_pets(overwrite, log_dir);
    case 5
        setup_leads_processing(data_dir, overwrite, process_unused_mris, cleanup);
        process_mris(overwrite, log_dir, segment_brainstem);
        process_pets(overwrite, log_dir);
    case 6
        fprintf('Ok\n');
    otherwise
        fprintf('INVALID INPUT. Try again :(\n');
end

% End the program
say_bye();
toc;
clear;
