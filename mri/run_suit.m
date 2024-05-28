function outfiles = run_suit(mri_dir, fid, overwrite)
    % Run the SUIT pipeline on a single MRI scan.
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % Setup parameters
    SUIT_VERSION = '3p7';
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);
    suit_dir = fullfile(mri_dir, 'suit');

    % Check if SUIT files already exist
    suitfs.t1 = fullfile(suit_dir, 't1.nii');  % The nu.nii in LPI orientation
    suitfs.gray = fullfile(suit_dir, 't1_seg1.nii');  % The cerebellar gray matter probability map in native MRI space
    suitfs.white = fullfile(suit_dir, 't1_seg2.nii');  % The cerebellar white matter probability map in native MRI space
    suitfs.isolation = fullfile(suit_dir, 'c_t1_pcereb.nii');  % The cerebellar mask in native MRI space
    suitfs.affineTr = fullfile(suit_dir, 'Affine_t1_seg1.mat');  % The affine transformation matrix
    suitfs.flowfield = fullfile(suit_dir, 'u_a_t1_seg1.nii');  % The flow field
    suitfs.iwAtlas = fullfile(suit_dir, 'iw_atl-Anatom_space-SUIT_dseg_u_a_t1_seg1.nii');  % The SUIT atlas in native MRI space
    suitfs.wgray = fullfile(suit_dir, 'wdt1_seg1.nii');  % The modulated cerebellar gray matter probability map in SUIT space
    suitfs.suit_vols = fullfile(suit_dir, append(scan_tag, '_cbl-suit-vols.csv')); % The SUIT ROI volumes
    outfiles.suit_atlas = fullfile(mri_dir, append(scan_tag, '_cbl-suit.nii')); % The SUIT atlas in RAS orientation

    if all(structfun(@isfile, outfiles)) && all(structfun(@isfile, suitfs)) && ~overwrite
        log_append(fid, '- SUIT files already exist, will not rerun');
        return
    else
        log_append(fid, '- Running SUIT');
    end

    % Create the SUIT directory
    suit_dir_full = fullfile(mri_dir, append('suit_', SUIT_VERSION));
    if ~isfolder(suit_dir_full)
        log_append(fid, '  * Creating the SUIT directory');
        mkdir(suit_dir_full);
        cmd = sprintf('ln -s %s %s', suit_dir_full, suit_dir);
        system(cmd);
    end

    % Convert the T1 from RAS to LPI orientation
    log_append(fid, '  * Converting RAS to LPI orientation');
    mri_convert = 'mri_convert --out_orientation LPI';
    fsfs = get_freesurfer_files(mri_dir);
    infile = fsfs.nu;
    outfile = suitfs.t1;
    if overwrite || ~isfile(outfile)
        cmd = sprintf('%s %s %s', mri_convert, infile, outfile);
        msg = sprintf( ...
            '    - %s ->\n                %s', ...
            basename(infile), ...
            basename(outfile) ...
        );
        log_append(fid, msg);
        run_system(cmd, 1, false, true);
    end

    % Segment the cerebellum
    log_append(fid, '  * Segmenting the cerebellum');
    if overwrite || ~all(isfile({suitfs.gray, suitfs.white, suitfs.isolation}));
        suit_isolate_seg({suitfs.t1}, 'maskp', 0.3);
    end

    % Normalize to the SUIT template
    log_append(fid, '  * Estimating transform from native MRI to SUIT space');
    if overwrite || ~all(isfile({suitfs.affineTr, suitfs.flowfield}))
        clear job;
        job.subjND(1).gray = {suitfs.gray};
        job.subjND(1).white = {suitfs.white};
        job.subjND(1).isolation = {suitfs.isolation};
        suit_normalize_dartel(job);
    end

    % Inverse warp the SUIT template to native MRI space
    log_append(fid, '  * Inverse warping the SUIT template to native MRI space');
    if overwrite || ~isfile(suitfs.iwAtlas)
        clear job;
        job.Affine = {suitfs.affineTr};
        job.flowfield = {suitfs.flowfield};
        job.ref = {suitfs.t1};
        suit_reslice_dartel_inv(job);
    end

    % Save a CSV file of the SUIT ROI volumes
    if overwrite || ~isfile(suitfs.suit_vols)
        suit_vols = suit_vol(suitfs.iwAtlas, 'Atlas');
        writetable(struct2table(suit_vols), suitfs.suit_vols);
    end

    % Convert SUIT files in subject space back to RAS orientation
    log_append(fid, '  * Converting back to the original orientation');
    if overwrite || ~isfile(outfiles.suit_atlas)
        % Change the datatype from int8 to int6 so mri_convert can handle it
        nii_change_datatype(suitfs.iwAtlas, spm_type('int16'));

        % Convert LPI to RAS orientation
        mri_convert = 'mri_convert --out_orientation RAS -rt nearest';
        infile = add_presuf(suitfs.iwAtlas, 'p');
        outfile = outfiles.suit_atlas;
        cmd = sprintf('%s %s %s', mri_convert, infile, outfile);
        msg = sprintf( ...
            '    - %s ->\n                %s', ...
            basename(infile), ...
            basename(outfile) ...
        );
        log_append(fid, msg);
        run_system(cmd, 1, false, true);

        % Reslice to match the original T1
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(fsfs.nu);
        matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(outfiles.suit_atlas);
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        spm_jobman('run',matlabbatch); clear matlabbatch;
        infile = add_presuf(outfiles.suit_atlas, 'r');
        outfile = outfiles.suit_atlas;
        movefile(infile, outfile);
    end

    % Warp cerebellar gray matter mask to SUIT space in preparation for VBM
    % analysis. GM probability is modulated based on the Jacobian
    log_append(fid, '  * Warping cerebellar GM probability map from native MRI to SUIT space');
    if overwrite || ~isfile(suitfs.wgray)
        clear job;
        job.subj.affineTr = {suitfs.affineTr};
        job.subj.flowfield = {suitfs.flowfield};
        job.subj.resample = {suitfs.gray};
        job.subj.jactransf = 1;
        job.subj.mask = {suitfs.isolation};
        suit_reslice_dartel(job);
    end
end
