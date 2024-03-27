# leads_processing
LEADS PET processing codebase

## Major changes in the 2024 update.

### Corrections
- Corrected an error in which nifti images saved with SPM's imcalc were
cast to an inappropriate data type (int16 for continuous values). These
images are now saved as float32
- Corrected an error in the algorithm to calculate a global amyloid SUVR
and corresponding Centiloid value using the ADNI longitudinal reference
region. Our original code calculated a reference region mean by taking
a weighted average across all voxels in the eroded WM, whole cerebellum,
and brainstem (combined into a single mask). But ADNI calculates a mean
within each of these 3 regions separately, then takes a straight average
across them, unweighted by volume.

### Organization and file naming
- LEADS project now located at `/mnt/coredata/processing/leads`
- Reorganized LEADS directory structure and file names in the following
ways:
    - Raw data is stored in `/mnt/coredata/processing/leads/data/raw`
    - Processed data is stored in `/mnt/coredata/processing/leads/data/processed`
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

### Processing
- All PET data now processed at 6mm resolution
- Dicoms now automatically converted to nifti using `dcm2niix`
- All reference regions are saved as masks in MRI native space at
`/mnt/coredata/processing/leads/data/processed/<subj>/MRI-T1_<YYYY-MM-DD>`
- The nu.nii file nonlinearly warped to template space is now saved as
`w<subj>_MRI-T1_<YYYY-MM-DD>_nu.nii`
