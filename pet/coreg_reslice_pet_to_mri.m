function outfile = coreg_reslice_pet_to_mri(infile, nuf, fid, overwrite, prefix)
    % Coregister and reslice PET scans to the MRI in the same directory.
    %
    % Parameters
    % ----------
    % infile : char|string
    %     Path to the PET image to coregister and reslice to MRI.
    % nuf : char|string
    %     Path to the MRI to coregister to. If empty, the MRI is looked
    %     for in <pet_dir>/mri
    % fid : int, optional
    %     File identifier for logging. Default is 1 for stdout
    % overwrite : logical
    %     If true, overwrite existing files. Default is false
    % prefix : char|string, optional
    %     Prefix to append to the output filenames. Default is 'r'
    %
    % Returns
    % -------
    % outfile : char|string
    %   Path to the coregistered and resliced PET image.
    % ------------------------------------------------------------------
    arguments
        infile {mustBeFile}
        nuf {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        prefix = 'r'
    end

    % Find the nu.nii
    if isempty(nuf)
        mri_dir = fullfile(dirname(infile), 'mri');
        mri_tag = get_scan_tag(mri_dir);
        nuf = fullfile(mri_dir, append(mri_tag, '_nu.nii'));
        mustBeFile(nuf);
    end

    % Format paths
    infile = abspath(infile);
    nuf = abspath(nuf);

    % Check existence of output files
    outfile = add_presuf(infile, prefix);
    if isfile(outfile) && ~overwrite
        log_append(fid, '- Coregistered PET image exists, will not overwrite');
        return
    else
        log_append(fid, '- Coregistering and reslicing PET to native MRI space:');
        if ~isempty(prefix)
            msg = sprintf( ...
                '  * %s ->\n              %s', ...
                basename(infile), ...
                basename(outfile) ...
            );
            log_append(fid, msg);
        else
            log_append(fid, sprintf('  * %s', basename(infile)));
        end
    end

    % Coregister and reslice PET to MRI
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {nuf};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {infile};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
end
