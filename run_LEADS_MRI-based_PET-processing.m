% Run the LEADS MRI-based PET processing pipeline with user input on the
% tasks to be completed. Full pipeline includes:
%
% Setup:
% - DICOM conversion to NIfTI
% - MRI and PET scans moved from newdata to raw directories
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
% - PET SUVR means calculcated for all FreeSurfer regions in the
%   aparc+aseg.nii + the brainstem sublabels file


%   always the mean of the cerebellar cortex
%   coregistered to the subject's MRI and then to the SPM12 OldNorm

% nu.nii MRI images to each subject's baseline MRI,
%  to  of
% all  SPM12  and
% SPM12, and PET processing to native-space MRI and MNI space.
%
%                                          Written by D. Schonhaut, 2024
% ----------------------------------------------------------------------

% Add paths to the scripts we need to access
restoredefaultpath;
addpath("/home/mac/ycobigo/neuroimaging/neuroimaging_CentOS7/SPM/spm12");
addpath("/home/mac/ycobigo/neuroimaging/neuroimaging_CentOS7/SPM/spm12/toolbox");
addpath(genpath(dirname(mfilename('fullpath'))));

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

% Ask the user what they want to do
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
say_bye();
clear;
