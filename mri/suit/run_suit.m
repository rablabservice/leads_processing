function outfiles = run_suit(mri_dir, overwrite)
    % Run the SUIT pipeline on a single MRI scan.
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
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
    suitfs.wgray = fullfile(suit_dir, 'wt1_seg1.nii');  % The modulated cerebellar gray matter probability map in SUIT space
    outfiles.suit_vols = fullfile(mri_dir, append(scan_tag, '_cbl-suit-vols.csv')); % The SUIT ROI volumes
    outfiles.suit_atlas = fullfile(mri_dir, append(scan_tag, '_cbl-suit.nii')); % The SUIT atlas in RAS orientation

    if all(structfun(@isfile, suitfs)) && ~overwrite
        fprintf('- SUIT files already exist, will not rerun\n');
        return
    else
        fprintf('- Running SUIT\n');
    end

    % Create the SUIT directory
    suit_dir_full = fullfile(mri_dir, append('suit_', SUIT_VERSION));
    if ~isfolder(suit_dir_full)
        fprintf('  * Creating the SUIT directory\n')
        mkdir(suit_dir_full);
        cmd = sprintf('ln -s %s %s', suit_dir_full, suit_dir);
        run_system_cmd(cmd, true, false);
    end

    % Convert the T1 from RAS to LPI orientation
    fprintf('  * Converting RAS to LPI orientation\n')
    mri_convert = 'mri_convert --out_orientation LPI';
    fsfs = get_freesurfer_files(mri_dir);
    infile = fsfs.nu;
    outfile = suitfs.t1;
    if overwrite || ~isfile(outfile)
        cmd = sprintf('%s %s %s', mri_convert, infile, outfile);
        fprintf('    - %s -> %s\n', basename(infile), basename(outfile));
        run_system_cmd(cmd, true, false);
    end

    % Segment the cerebellum
    fprintf('  * Segmenting the cerebellum\n');
    if overwrite || ~all(isfile({suitfs.gray, suitfs.white, suitfs.isolation}))
        suit_isolate_seg({suitfs.t1}, 'maskp', 0.3);
    end

    % Normalize to the SUIT template
    fprintf('  * Estimating transform from native MRI to SUIT space\n');
    if overwrite || ~all(isfile({suitfs.affineTr, suitfs.flowfield}))
        clear job;
        job.subjND(1).gray = {suitfs.gray};
        job.subjND(1).white = {suitfs.white};
        job.subjND(1).isolation = {suitfs.isolation};
        suit_normalize_dartel(job);
    end

    % Inverse warp the SUIT template to native MRI space
    fprintf('  * Inverse warping the SUIT template to native MRI space\n');
    if overwrite || ~isfile(suitfs.iwAtlas)
        clear job;
        job.Affine = {suitfs.affineTr};
        job.flowfield = {suitfs.flowfield};
        job.ref = {suitfs.t1};
        suit_reslice_dartel_inv(job);
    end

    % Save a CSV file of the SUIT ROI volumes
    if overwrite || ~isfile(outfiles.suit_vols)
        fprintf('  * Retrieving cerebellar lobule volumes\n');
        suit_vols = suit_vol(suitfs.iwAtlas, 'Atlas');
        writetable(struct2table(suit_vols), outfiles.suit_vols);
        fprintf('    - Saved %s\n', basename(outfiles.suit_vols));
    end

    % Convert SUIT files in subject space back to RAS orientation
    fprintf('  * Converting back to the original orientation\n')
    if overwrite || ~isfile(outfiles.suit_atlas)
        % Change the datatype from int8 to int6 so mri_convert can handle it
        nii_change_datatype(suitfs.iwAtlas, spm_type('int16'));

        % Convert LPI to RAS orientation
        mri_convert = 'mri_convert --out_orientation RAS -rt nearest';
        infile = add_presuf(suitfs.iwAtlas, 'p');
        outfile = outfiles.suit_atlas;
        cmd = sprintf('%s %s %s', mri_convert, infile, outfile);
        fprintf('    - %s -> %s\n', basename(infile), basename(outfile));
        run_system_cmd(cmd, true, false);
        delete(infile);

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
    fprintf('  * Warping cerebellar GM probability map from native MRI to SUIT space\n');
    if overwrite || ~isfile(suitfs.iwAtlas)
        clear job;
        job.subj.affineTr = {suitfs.affineTr};
        job.subj.flowfield = {suitfs.flowfield};
        job.subj.resample = {suitfs.gray};
        job.subj.jactransf = 1;
        job.subj.mask = {suitfs.isolation};
        suit_reslice_dartel(job);
    end
end
