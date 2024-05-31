# LEADS MRI-based PET processing
Data processing pipeline for the Longitudinal Early-onset Alzheimer's Disease Study (LEADS).

## Primer
This repository contains the new codebase to process LEADS MRI and PET scans using an
MRI-based processing pipeline that mimics the ADNI processing pipeline. The new codebase
was written by D. Schonhaut in Spring 2024 as a full rewrite of the original LEADS processing
code by L. Iaccarino in 2018.

The new codebase is designed to be faster, more readable, and more modular than the
original codebase. With minimal changes, the new code should also be fully extesible to other
projects that collect PET and MRI data and want to adopt an MRI-based PET processing pipeline
with a similar directory structure and file naming conventions.

## Corrections to the original LEADS processing pipeline
The points below encompass, to the best of our knowledge, all *substantive* differences
between the new LEADS pipeline and the older processing code it replaces. Non-substantive
differences refer to any changes that do not affect the underlying values of processed data
(for example, changes to file naming convenctions and directory organization, or changes to
the code itself that don't affect the end result for files saved).
1. The original codebase used SPM ImCalc to perform arithmetic operations on PET images
   (e.g. create an SUVR by dividing PET voxels by a scalar value). However, the ImCalc code
   was programmed to cast output values to an inappropriate data type (int16 for continuous
   values) and reduced the decimal precision vs. the input values by 4-fold. In the new
   codebase, SPM ImCalc is no longer used (working with Matlab arrays is faster) but all
   NIfTI images that hold continuous values are stored at float32 precision.
1. Corrected an error in the algorithm to calculate a global amyloid SUVR and
   corresponding Centiloid value using the ADNI longitudinal reference region. To obtain
   the reference region scaling factor, our original code calculated a weighted mean across
   all voxels in the eroded WM, whole cerebellum, and brainstem after combining these
   regions into a single mask. In fact, ADNI calculates the longitudinal reference region
   scaling factor by calculating the mean PET signal within the eroded WM, whole
   cerebellum, and brainste separately, then taking a straight average across these
   regions, unweighted by their respective volumes.
1. To create the inferior cerebellar gray matter mask that is used as a reference region
   for tau-PET, the cerebellar SUIT atlas is warped into the subject's native MRI space. In
   the original LEADS code, warped SUIT atlas files were resliced with trilinear
   interpolation and then re-cast as integer values. This caused edge effects with non-
   sensical values at the borders between different cerebellar subregions (e.g. a left
   subregion label appearing in R hemisphere). The new code fixes this problem by using
   nearest-neighbor interpolation to reslice the SUIT atlas to native MRI space.
1. Also regarding to the inferior cerebrellar reference region, we now use the
   [SUIT toolbox](https://www.diedrichsenlab.org/imaging/suit.htm) in MATLAB to estimate
   a nonlinear transformation between each subject's cerebellum in native space and the
   'SUIT space' template cerebellum. Previously we had used the iy_ file defined during
   SPM segmentation to reverse normalize an MNI space version of the SUIT template.
1. The regularization parameters used for SPM segmentation in the original LEADS codebase
   differed slightly from the default regularization parameters used by SPM12.
   (original LEADS: matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.025 0.1];)
   (default params: matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];).
   We could not find a rationale for this change and did not observe a noticeable
   difference in segmentation results when testing segmentation with the altered versus
   default parameters. For simplicity, the new codebase reverts back to the SPM12 default.
1. The original LEADS codebase saved PET images that were affine transformed to MNI space
   at 1mm^3 resolution, while PET images that were nonlinearly warped to MNI space were
   saved at 1.5mm^3. In the new LEADS codebase, affine and warped MNI space images are both
   saved at 1.5mm^3.

## Structure and implementation of the new codebase
The new codebase is implemented in Matlab (version 2023b) and Python (version 3.11), with
Matlab serving as the primary interface for end users. The user should not need to specify
paths to any programs aside from Matlab, as all paths are specified by the code to ensure
that the same versions of required software are used consistently.

A few additional notes about the new codebase:
- On the PETcore VM, the codebase lives in `/mnt/coredata/processing/leads/code` and is
  connected to its [remote repository](https://github.com/dschonhaut/leads_processing) on GitHub.
- The new code is 10-20x faster than the original codebase, mostly due to the use of multithreading
  throughout the processing pipeline. FreeSurfer still takes ~10 hr to run per scan, although up to
  16 MRIs can be processed in parallel on the PETcore VM.
- A single Matlab script, [`run_leads_mri_based_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/run/run_leads_mri_based_processing.m), allows
  the user to run all parts of the processing pipeline while interacting with the program
  through simple, question and answer based prompts.
- The processing pipeline is divided into three modules, one of which must be specified by
  the user with each call to [`run_leads_mri_based_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/run/run_leads_mri_based_processing.m).

## Running the processing pipeline
From an end user's perspective, new MRI and PET data will be downloaded from LONI and manually
uploaded as one or more zip files to `/mnt/coredata/processing/leads/data/newdata`. The
`run_leads_mri_based_processing.m` program then needs to be run four times in the following
order:
1. Run the **Setup Module** (Option 1, "Setup scans for processing") to move scans from
   `.../leads/data/newdata` to `.../leads/data/raw` and schedule new MRIs for processing
   (**_note_**: scans are selected but not yet processed). This module takes 2-3 min to
   run, and at the end the user is informed of how many scans have been scheduled for
   processing.
   1. Optionally, the user can now run `run_leads_mri_based_processing.m` and choose Option 2
      ("View scans that are scheduled for processing") to see a full list of scheduled scans.
1. Run the **MRI Module** (Option 3, "Process scheduled MRIs") to process new MRIs through
   FreeSurfer and SPM-based pipelines.
1. Rerun the **Setup Module** to schedule new PET scans for processing. This step is necessary
   because PET scans will not be scheduled until the MRI that they will be coregistered to has
   been processed.
   1. As before, running `run_leads_mri_based_processing.m` and choosing Option 2 ("View scans
      that are scheduled for processing") will print a full list of scheduled PET scans.
1. Run the **PET Module** (Option 4, "Process scheduled PET scans") to process new PET scans.

## Data processing modules
Here is a somewhat more detailed account of what happens during each processing module.

### Setup
* Main script: [`setup_leads_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/setup/setup_leads_processing.m)
* Overview: Newly uploaded data are prepared for subsequent processing
    1. Zip files in `.../leads/data/newdata` are unzipped
    1. All DICOMs in `.../leads/data/newdata` are converted to NIfTI using
    [`dcm2niix`](https://github.com/rordenlab/dcm2niix) (version v1.0.20240202).
    1. MRI and PET scan directories are moved from `.../leads/data/newdata` to `.../leads/data/raw`
    if they don't already exist in `.../leads/data/raw`
    1. All scan directories in `.../leads/data/raw` are identified by a recursive search for
    *.nii files
    1. Filenames are parsed to resolve:
    -  Subject ID
    -  Scan type (MRI or PET; and for PET, which tracer)
    -  Scan acquisition date
    -  LONI Image ID
    -  For PET scans, the spatial resolution
    1. Each PET scan is matched to the closest MRI scan
    1. Orphan MRIs (those not recently acquired and not closest to any PET
    scan) are identified, and under default settings will not be scheduled
    for processing
    1. PET scans are audited for potential issues that would prevent processing.
    These include:
    - Inability to parse any of the scan information above
    - PET not at the expected resolution (6mm by default)
    - PET and MRI scans too far apart (>365 days, by default)
    - Same MRI would be used for multiple PET scans of the same tracer (e.g.
        two FTP timepoints)
    1. A search is conducted to identify which MRIs in `.../leads/data/raw` have been
    processed by FreeSurfer, and which have completed post-FreeSurfer,
    SPM-based processing, respectively
    1. A search is conducted to identify which PET scans in `.../leads/data/raw` have been
    fully processed
    1. MRIs are scheduled for processing if they have not been fully processed,
    are not orphans, and (for follow-up MRIs) their corresponding baseline
    MRI has been processed
    1. PET scans are scheduled for processing if they have not been fully
    processed, are not flagged with an issue, and their corresponding MRI
    has been fully processed
    1. The user is informed of how many MRI and PET scans are in `.../leads/data/raw`,
    and how many of these scans have completed processing, been flagged with one or
    more issues, and been scheduled for processing, respectively
    1. Two CSV files containing scan-level information are saved to
    'scans_to_process': one for MRIs and one for PET scans. These files are
    read by downstream scripts in the processing pipeline. Manually editing
    these files will affect which scans the scripts attempt to process,
    though could also introduce errors, so edit the CSVs with caution and at
    your own risk.
    1. For scans that are scheduled to be processed, new directories are created in
    `.../leads/data/processed` to store the processed data. Symbolic links are made in
    these new scan directories that point back to their corresponding directories
    in `.../leads/data/raw`. A symbolic link in each PET scan directory is also created that
    points to the directory of the MRI it will be coregistered to.

### MRI
* Main script: [`process_mris.m`](https://github.com/dschonhaut/leads_processing/blob/main/mri/process_mris.m)
* Overview: MRI processing is divided into two submodules:
    1. T1 MRIs are processed through FreeSurfer [`recon-all`](https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all)
       and [`segmentBS.sh`](https://surfer.nmr.mgh.harvard.edu/fswiki/BrainstemSubstructures), which
       segments the brainstem into subregions. This is the most time consuming part of the processing
       pipeline, taking ~10 hours to run on a single MRI, and currently the PETcore VM is equipped
       to process only up to 16 scans in parallel. In contrast, the second most time consuming part
       of processing is SPM segmentation, which takes about 12 minutes per scan.
    1. Post-FreeSurfer processing is carried out in [`SPM12`](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).
       The steps in order are:
       1. FreeSurfer files are converted from .mgz to .nii and copied from the FreeSurfer directory
          into the processed MRI directory
       1. Image origins of the nu.nii, aparc+aseg.nii, and brainstem_sublabels.nii files are reset to the
          center-of-mass of the nu.nii. Rigid body coregistration parameters are then estimated to bring
          the nu.nii as close as possible to the SPM12 OldNorm T1.nii template in MNI space, and these
          parameters are applied to the nu.nii, aparc+aseg.nii, and brainstem_sublabels.nii images.
       1. The nu.nii is coregistered to the baseline nu.nii (earliest processed MRI for a given subject),
          and the same transform is applied to the aparc+aseg.nii and brainstem_sublabels.nii images.
          This step is skipped for baseline MRIs.
       1. The [Matlab SUIT toolbox](https://www.diedrichsenlab.org/imaging/suit.htm) is run to estimate
          a nonlinear fit between the native space nu.nii and the cerebellar SUIT template.
       1. Reference region and target ROI masks are saved in native MRI space.
       1. The nu.nii is segmented into tissue probability maps, and deformation fields that map
          the nonlinear transforms between native MRI and MNI space are saved (y_* and iy_* files).
       1. The nu.nii is warped to MNI space using the y_* deformation field
       1. The full (12-degree) affine transform between the native space nu.nii and the OldNorm
          T1.nii template in MNI space is estimated and saved.
       1. The affine transform from the previous step is applied to the native space nu.nii

### PET
* Main script: [`process_pets.m`](https://github.com/dschonhaut/leads_processing/blob/main/pet/process_pets.m)
* Overview: PET processing must be performed after the MRI that PET will be coregistered to has
    been processed, and the [Setup module](#setup) will not schedule PET scans to be processed until this
    has happened. PET processing includes the following steps, in order:
    1. The "raw" PET scan is copied from `.../leads/data/raw` to `.../leads/data/processed`
       - **_Note:_** In LEADS we download already preprocessed PET data, for which reconstructed PET frames have
         already been realigned and averaged into a single frame spanning a tracer-specific time window of
         interest. Preprocessed scans have also been resliced to the same voxel size, intensity normalized against
         mean signal in the whole cerebellum (not FreeSurfer-based), and smoothed to approximately the same resolution.
         This code will therefore need to be modified before it can be used to process raw PET data (reconstructed frames)
         rather than preprocessed PET data.
    1. The PET scan is coregistered and resliced to the nu.nii (at 1mm^3 voxel sizes)
    1. Reference region means are calculated, and voxelwise PET SUVR images are saved.
       - Reference regions in the new pipeline are not hard-coded, but rather depend on a
         [`ref_regions.csv`](https://github.com/dschonhaut/leads_processing/blob/main/config/ref_regions.csv)
         file that in the `.../leads/code/config` directory. To process a new PET tracer or reference region,
         the user should need only to add a line to this file and -- if the reference region is not already
         saved in the processed MRI directory -- add a new script to save the desired reference region
         and call it from [`save_roi_masks.m`](https://github.com/dschonhaut/leads_processing/blob/main/mri/save_roi_masks.m).
         If the PET tracer is not recognized, it may need to be added to a different config file,
         [`scan_types_and_tracers.csv`](https://github.com/dschonhaut/leads_processing/blob/main/config/scan_types_and_tracers.csv).
       - Following ADNI processing, the default reference regions saved in LEADS are:
         * Amyloid-PET
           - Whole cerebellum (cross-sectional reference region)
           - Composite white matter, consisting of the unweighted average PET signal within the brainstem,
             eroded subcortical white matter, and whole cerebellum (longitudinal reference region)
         * Tau-PET
           - Inferior cerebellar gray matter (cross-sectional reference region)
           - Eroded subcortical white matter (longitudinal reference region)
         * FDG-PET
           - Pons
    1. Mean PET SUVRs and ROI volumes are extracted from FreeSurfer regions in native MRI space.
      - The `fsroi_list_<TRACER>.csv` files in
        [`../leads/code/config`](https://github.com/dschonhaut/leads_processing/tree/main/config)
        define which ROIs are saved for each PET tracer
    1. For amyloid-PET only, the mean SUVR within the ADNI cortical summary region is calculated
       and converted to Centiloids using the appropriate equation from
       [Royse et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33971965)
    1. PET SUVRs are warped from native MRI space to MNI space using the y_* deformation field
       obtained during segmentation of the nu.nii.
    1. PET SUVRs are linearly transformed from native MRI space to MNI space using the affine transform
       estimated on the nu.nii.
