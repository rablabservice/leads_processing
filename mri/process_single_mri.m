function outfiles = process_single_mri( ...
    mri_dir, ...
    raw_mrif, ...
    overwrite, ...
    segment_brainstem, ...
    process_freesurfer, ...
    process_post_freesurfer ...
)
    % Run a single MRI scan through all the processing steps.
    %
    % Overview
    % --------
    % 1.  Run recon-all on the raw T1 MRI
    % 2.  Optionally segment the brainstem using segmentBS.sh
    % 3.  Copy and convert FreeSurfer files to NIfTI
    % 4.  Reset FreeSurfer file origins to the nu.nii center-of-mass,
    %     then coregister nu.nii to MNI space
    % 5.  Coregister nu.nii to the baseline nu.nii (earliest MRI for a
    %     subject). This step is skipped if nu.nii is at baseline.
    %     FreeSurfer files are now in native MRI space
    % 6.  Run SUIT and inverse warp cerebellar subregions to native MRI
    %     space
    % 7.  Save reference region and target ROI masks in native MRI space
    % 8.  Segment nu.nii and save deformation fields between native MRI
    %     and MNI space
    % 9.  Warp nu.nii to MNI space
    % 10. Use nu.nii to calculate and save the full affine transform
    %     from native MRI to MNI space
    % 11. Apply affine transformation to the native space nu.nii
    %
    % Usage
    % -----
    % outfiles = process_single_mri(mri_dir)
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %     Path to the processed MRI directory
    % raw_mrif : char or str array
    %     Full path to the raw MRI nifti. If not passed, this is assumed
    %     to be the first .nii file in mri_dir/raw
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is false
    % segment_brainstem : logical, optional
    %     If true, segment the brainstem using segmentBS.sh. Default is
    %     true
    % process_freesurfer : logical, optional
    %     If true, run recon-all on the raw MRI. Default is true
    % process_post_freesurfer : logical, optional
    %     If true, run post-FreeSurfer processing. Default is true
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        raw_mrif {mustBeText} = ''
        overwrite logical = false
        segment_brainstem logical = true
        process_freesurfer logical = true
        process_post_freesurfer logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % If processing is already complete and overwrite is false, get
    % the struct of processed MRI files and return
    if processed_mri_files_exist(mri_dir) && ~overwrite
        fprintf('%s processing already complete, returning output files\n', scan_tag)
        outfiles = get_processed_mri_files(mri_dir);
        return
    end

    % If raw_mrif is not passed, find in it mri_dir/raw
    if isempty(raw_mrif)
        files = dir(fullfile(mri_dir, 'raw', '*.nii'));
        if isempty(files)
            error('No raw MRI files found in %s', fullfile(mri_dir, 'raw'));
        elseif length(files) > 1
            error('Multiple raw MRI files found in %s', fullfile(mri_dir, 'raw'));
        end
        raw_mrif = abspath(fullfile({files.folder}, {files.name}));
    end

    % Start the log file
    fid = log_start(mri_dir);

    % Print the module header
    title = 'MRI PROCESSING MODULE';
    subtitle = append('SCAN : ', scan_tag);
    print_header(title, subtitle, fid);

    % Run FreeSurfer
    if process_freesurfer
        process_mri_freesurfer(raw_mrif, mri_dir, fid, overwrite, segment_brainstem);
    end

    % Run post-FreeSurfer processing
    if process_post_freesurfer
        if freesurfer_files_exist(mri_dir)
            outfiles = process_mri_post_freesurfer(mri_dir, fid, overwrite);
        else
            log_append(fid, '- Cannot complete post-FreeSurfer processing until all FreeSurfer files exist');
        end
    end

    % Print the module footer
    print_footer('MRI processing module complete', fid);

    % Close the log file
    log_close(fid);
end
