function outfiles = process_single_mri( ...
    mri_dir, ...
    overwrite, ...
    raw_mrif, ...
    segment_brainstem, ...
    process_freesurfer, ...
    process_post_freesurfer, ...
    run_qc ...
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
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is false
    % raw_mrif : char or str array
    %     Full path to the raw MRI nifti. If not passed, this is assumed
    %     to be the first .nii file in mri_dir/raw
    % segment_brainstem : logical, optional
    %     If true, segment the brainstem using segmentBS.sh. Default is
    %     true
    % process_freesurfer : logical, optional
    %     If true, run recon-all on the raw MRI. Default is true
    % process_post_freesurfer : logical, optional
    %     If true, run post-FreeSurfer processing. Default is true
    % run_qc : logical, optional
    %     If true, run the QC script. Default is true
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        raw_mrif {mustBeText} = ''
        segment_brainstem logical = true
        process_freesurfer logical = true
        process_post_freesurfer logical = true
        run_qc logical = true
    end

    % ------------------------------------------------------------------
    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % ------------------------------------------------------------------
    % If processing is already complete and overwrite is false, get
    % the struct of processed MRI files and return
    if processed_mri_files_exist(mri_dir) && ~overwrite
        fprintf('%s processing already complete, returning output files\n', scan_tag)
        outfiles = get_processed_mri_files(mri_dir);
        return
    end

    % ------------------------------------------------------------------
    % If raw_mrif is not passed, find in it mri_dir/raw
    if isempty(raw_mrif)
        files = dir(fullfile(mri_dir, 'raw', '*.nii'));
        if isempty(files)
            error('No raw MRI files found in %s', fullfile(mri_dir, 'raw'));
        elseif length(files) > 1
            error('Multiple raw MRI files found in %s', fullfile(mri_dir, 'raw'));
        end
        raw_mrif = abspath(fullfile({files.folder}, {files.name}));
        raw_mrif = raw_mrif{1};
    end

    % ------------------------------------------------------------------
    % Start the log file
    fid = log_start(mri_dir);

    try
        % --------------------------------------------------------------
        % Print the module header
        log_append(fid, '', 0, 0);
        log_append(fid, 'START MRI PROCESSING MODULE', 0, 0);
        log_append(fid, '---------------------------', 0, 0);
        log_append(fid, '', 0, 0);

        % Print path to the current code file
        log_append(fid, 'Code file:', 0, 0);
        log_append(fid, sprintf('%s.m\n', mfilename('fullpath')), 0, 0);

        % Print the input parameters
        log_append(fid, 'Input parameters:', 0, 0);
        log_append(fid, sprintf('mri_dir = %s', mri_dir), 0, 0);
        log_append(fid, sprintf('overwrite = %d', overwrite), 0, 0);
        if isempty(raw_mrif)
            log_append(fid, 'raw_mrif = ''''', 0, 0);
        else
            log_append(fid, sprintf('raw_mrif = %s', raw_mrif), 0, 0);
        end
        log_append(fid, sprintf('segment_brainstem = %d', segment_brainstem), 0, 0);
        log_append(fid, sprintf('process_freesurfer = %d', process_freesurfer), 0, 0);
        log_append(fid, sprintf('process_post_freesurfer = %d', process_post_freesurfer), 0, 0);
        log_append(fid, '', 0, 0);

        % --------------------------------------------------------------
        % Run FreeSurfer
        if process_freesurfer
            process_mri_freesurfer(raw_mrif, mri_dir, fid, overwrite, segment_brainstem);
        end

        % --------------------------------------------------------------
        % Run post-FreeSurfer processing
        if process_post_freesurfer
            if freesurfer_files_exist(mri_dir)
                outfiles = process_mri_post_freesurfer(mri_dir, fid, overwrite);
            else
                log_append(fid, '- Cannot complete post-FreeSurfer processing until all FreeSurfer files exist');
            end
        end

        % --------------------------------------------------------------
        % Save the QC image
        if run_qc
            log_append(fid, '- Generating QC image');
            python = '/home/mac/dschonhaut/mambaforge/envs/nipy311/bin/python';
            code_dir = fileparts(fileparts(mfilename('fullpath')));
            qc_script = fullfile(code_dir, 'qc', 'leadsqc.py');
            cmd = sprintf('%s %s %s', python, qc_script, mri_dir);
            system(cmd);
        end

        % --------------------------------------------------------------
        % Identify and log any missing processed files
        log_append(fid, '- Checking for expected output files');
        if processed_mri_files_exist(mri_dir)
            log_append(fid, '  * Found all expected output files');
        else
            outfiles_expected = cellvec(get_processed_mri_files(mri_dir));
            outfiles_missing = outfiles_expected(~cellfun(@isfile, outfiles_expected));
            n_missing = length(outfiles_missing);
            if n_missing == 1
                log_append(fid, '  * WARNING: Missing 1 expected output');
            else
                log_append( ...
                    fid, ...
                    sprintf('  * WARNING: Missing %d expected outputs', n_missing) ...
                );
            end
            cellfun(@(x) log_append(fid, sprintf('    - %s', x)), outfiles_missing);
        end

        % --------------------------------------------------------------
        % Print the module footer
        log_append(fid, '', 0, 0);
        log_append(fid, '-------------------------', 0, 0);
        log_append(fid, 'END MRI PROCESSING MODULE', 0, 0);
        log_append(fid, '', 0, 0);

        % --------------------------------------------------------------
        % Close the log file
        log_close(fid);
    catch ME
        % Print the error message
        log_append(fid, '!! ERROR !!');
        log_append(fid, getReport(ME, 'extended', 'hyperlinks', 'off'), 0, 0);

        % Close the log file
        log_append(fid, '\nClosing log file early due to error', 0, 0);
        log_close(fid);

        % Rethrow the error
        rethrow(ME);
    end
end
