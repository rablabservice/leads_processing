# LEADS MRI-based PET processing
Data processing pipeline for the Longitudinal Early-onset Alzheimer's Disease Study (LEADS).

## Primer
This repository contains the updated (version 2) codebase to process LEADS MRI and PET
scans under the MRI-based processing pipeline. The new LEADS codebase, introduced in Spring 2024, is a full rewrite
of the original MRI-based PET processing pipeline that was developed in 2018. Both codebases mimic ADNI processing, and aside from the corrections and changes noted in this document, the new codebase follows the methods used by the original pipeline. The new codebase is designed to be faster,
more readable, more modular, and more flexible than the original codebase. With minor changes it should be fully extensible to other projects that want to adopt an MRI-based PET processing pipeline using a similar directory structure and file naming conventions.

## Corrections
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

## Structure and implementation of the new codebase
- The new codebase is implemented in Matlab and Python, with Matlab as the primary
  interface. Aside from launching Matlab (the code was written in Matlab version 2023b), the
  user should not need to set paths to any other programs, as this is all done internally to
  ensure that all users are using the same versions of required software.
- On the PETcore VM, the codebase lives in **`/mnt/coredata/processing/leads/code`** and is
  connected to its remote repository on GitHub.
- The new codebase offers 10-20x performance improvements over the original codebase when
  run on the PETcore VM. This is mostly due to the use of multithreading throughout the
  processing pipeline, as opposed to only during FreeSurfer processing.
- A single Matlab script, [`run_leads_mri_based_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/run/run_leads_mri_based_processing.m), allows
  the user to run all parts of the processing pipeline while interacting with the program
  through simple, question and answer based prompts.
- The processing pipeline is divided into three modules, one of which must be specified by
  the user with each call to [`run_leads_mri_based_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/run/run_leads_mri_based_processing.m). The modules are:
    1. ### Setup
       * Main script: [`setup_leads_processing.m`](https://github.com/dschonhaut/leads_processing/blob/main/setup/setup_leads_processing.m)
       * Overview: Prepare newly downloaded data to be processed. There are three substages
         to this module:
         1. Newly downloaded files that the user copies to **`.../leads/data/newdata`** are unzipped,
            converted from DICOMs to NIfTI, and moved to **`.../leads/data/raw`**.
         1. All scans in **`.../leads/data/raw`** are parsed and compared against scans in
            **`.../leads/data/processed`** to determine which scans need to be processed. The results
            of this evaluation step are stored in two CSV files, *raw_MRI_index\*.csv* and
            *raw_PET_index\*.csv* in **`.../leads/metadata/scans_to_process`**. These files contain an
            entry for every MRI and PET scan in **`.../leads/data/raw`**, respectively, and the column
            "scheduled_for_processing" controls which scans the processing pipeline will attempt
            to process under the default settings.
         1. For scans that are scheduled to be processed, the code will create new directories in
            **`.../leads/data/processed`** to store the processed data. Symbolic links are made in
            these new, processed scan directories that point back to their corresponding directories
            in **`.../leads/data/raw`**. A symbolic link in each PET scan directory is also created that
            points to the directory of the closest MRI (the MRI scan that PET will be coregistered to).
    1. ### MRI processing
       * Main script: [`process_mris.m`](https://github.com/dschonhaut/leads_processing/blob/main/mri/process_mris.m)
       * Overview: MRI processing is divided into two submodules:
         1. T1 MRIs are processed through FreeSurfer [`recon-all`](https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all) and optionally [`segmentBS.sh`](https://surfer.nmr.mgh.harvard.edu/fswiki/BrainstemSubstructures), which
            segments the brainstem into subregions. This is the most time consuming part of the processing
            pipeline, taking ~10 hours to run on a single MRI, and currently the PETcore VM is equipped
            to process only up to 16 scans in parallel. In contrast, the second most time consuming part
            of processing is SPM segmentation, which takes about 12 minutes per scan.
         1. Post-FreeSurfer processing is carried out in [`SPM12`](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/). The steps in order are:
            1. Convert FreeSurfer files from .mgz to .nii, and copy them from the FreeSurfer directory
               to the processed MRI directory for a given scan.
            1. Reset origins of the nu.nii, aparc+aseg.nii, and brainstem_sublabels.nii images to the
               center-of-mass of the nu.nii. Then coregister the nu.nii to the OldNorm T1.nii template
               in MNI space, and apply the estimated transform to all three images.
            1. Coregister the nu.nii to the baseline nu.nii (earliest processed MRI for a given subject),
               and apply the estimated transform to the nu.nii, aparc+aseg.nii, and brainstem_sublabels.nii
               images. This step is skipped for baseline MRIs.
            1. Run the [Matlab SUIT toolbox](https://www.diedrichsenlab.org/imaging/suit.htm), which estimates
               a nonlinear fit between the cerebellar SUIT template and a native space MRI.
            1. Save reference region and target ROI masks in native MRI space.
            1. Segment the nu.nii and save deformation fields that map the nonlinear transform from native
               MRI to MNI space, and vice versa (y_* and iy_* files).
            1. Warp the nu.nii to MNI space using the y_* deformation field
            1. Calculate the full (12-degree) affine transform between the native space nu.nii and the OldNorm
               T1.nii template in MNI space, and save the affine transformation matrix.
            1. Apply the affine transformation to the native space nu.nii
    1. ### PET processing
       * Main script: [`process_pets.m`](https://github.com/dschonhaut/leads_processing/blob/main/pet/process_pets.m)
       * Overview: PET processing must be performed after the MRI that PET will be coregistered to has
         been processed, and the [Setup module](#setup) will not schedule PET scans to be processed until this
         is the case. PET processing includes the following steps, in order:
         1. Copy the PET scan from **`.../leads/data/raw`** to **`.../leads/data/processed`**
         1. Coregister and reslice PET to the nu.nii
         1. Calculate reference region means, and save voxelwise PET SUVR images. The reference regions
            used to make SUVR imagess are defined in [`ref_regions.csv`](https://github.com/dschonhaut/leads_processing/blob/main/config/ref_regions.csv).
            Following ADNI processing, default reference regions are:
            - Amyloid-PET
              * Whole cerebellum (cross-sectional reference region)
              * Composite white matter, consisting of the unweighted average PET signal within the brainstem,
                eroded subcortical white matter, and whole cerebellum (longitudinal reference region)
            - Tau-PET
              * Inferior cerebellar gray matter (cross-sectional reference region)
              * Eroded subcortical white matter (longitudinal reference region)
            - FDG-PET
              * Pons
         1. Extract mean PET SUVRs and ROI volumes within FreeSurfer regios in native MRI space.
            `fsroi_list_<TRACER>.csv` files in the [**`config`**](https://github.com/dschonhaut/leads_processing/tree/main/config) directory define which regions are saved for each PET tracer.
         1. For amyloid-PET only, calculate the mean SUVR within the cortical summary region
            and convert this value to Centiloids.
         1. Warp PET SUVRs to MNI space using the y_* deformation field from SPM12 segmentation
            of the nu.nii.
         1. Linearly transform PET SUVRs to MNI space using the affine transform estimated to
            scale the nu.nii to MNI space.

## Processing notes
- All PET data are now processed at 6mm resolution. The code checks for the resolution of
  PET scans in the downloaded filenames and alerts the user if the resolution is not 6mm
  or cannot be determined from the filename. Such scans are not, by default, processed
  without explicit direction from the user.
- Zip files downloaded from LONI are now automatically unzipped by the using a faster algorithm
  than the default `unzip` program in Linux.
- Dicoms are now automatically converted to nifti using [`dcm2niix`](https://github.com/rordenlab/dcm2niix) version v1.0.20240202.
-

### Data organization and file naming
1. The LEADS project now located at `/mnt/coredata/processing/leads`
1. A new Linux group, 'leads', has been created to manage write access to this project.
   All PETcore VM users retain the ability to read and copy out files within 'leads', but
   only group members will have write access to the data directories, and only I have
   write access to the codebase. This should offer greater security against files being
   accidentally modified or deleted by users who mean to copy it elsewhere.
1. The LEADS directory structure and file naming conventions have been reconsidered:
    - Raw data is stored in `/mnt/coredata/processing/leads/data/raw`
    - Processed data is stored in `/mnt/coredata/processing/leads/data/processed`
    - `MRI_T1` filenames and directories are being changed to `MRI-T1`
    - MRI-based mask files are being moved from PET directories to the
    processed `MRI-T1` directory, and PET directories will then symlink
    to the mask files that are used e.g. in making SUVR images
    - "w_affine" files are now simply prefixed with "a" to distinguish
    the linear affine transform from the nonlinear warping prefix "w",
    e.g. `a<subj>_MRI-T1_<YYYY-MM-DD>_nu.nii`
- Introduced a Python script (`LEADS_find_new_scans.py`) that identifies
scans that need to be processed by comparing unprocessed data that has
been uploaded to `/mnt/coredata/processing/leads/data/raw` against
processed data in `/mnt/coredata/processing/leads/data/processed`
    - This includes code that automatically parses the subject ID, scan
    date, and MRI sequence or tracer type from the file and directory
    names assigned by LONI and present at the time that the data are
    downloaded. This code makes use of a new spreadsheet,
    `/mnt/coredata/processing/leads/metadata/ssheets/scan_types_and_tracers.csv`,
    that maps LONI-assigned PET tracer names to whatever we want to
    call them in `/mnt/coredata/processing/leads/data/processed`
- LEADS code no longer requires subject IDs to begin with LDS*, but
instead identifies subjects by the directory name in which their data
is stored in `/mnt/coredata/processing/leads/data/raw`
- There are no more "Timepoint" directories in the processed file
hierarchy. Instead, each processed MRI and PET scan has its own
directory like `/mnt/coredata/processing/leads/data/processed/<subj>/<scan>`,
and PET scans use symbolic links to point to the MRIs they are
coregistered to.
- Hard-coded paths to template files at `/mnt/neuroimaging` have been
remapped to function outside of the old Singularity image
- Code to generate QC images has been moved out of the main processing
pipeline and will be implemented as a separate script
- Raw data is no longer removed at the end of processing
- Rewrote the LEADS extraction script and made several changes to the
output spreadsheets:
- The entire codebase has been refactored to trim length, provide
clarifying comments where needed, replace opaque variable names, reduce
hard-coding (especially in referencing dataframe columns by number or
substrings by position), create a more modular structure by making
better use of functions, and improve overall readability


