function outfiles = coreg_mri_to_baseline(infiles, baseline_nuf, prefix, overwrite)
    % Coregister MRI to baseline and ovewrite input image headers
    %
    % First looks up if infiles already correspond to the baseline
    % MRI, in which case nothing is done.
    %
    % Parameters
    % ----------
    % infiles : char/str or cell array of char/str
    %     Cell array of paths to the MRI files to coregister. The first
    %     element should be the MRI to coregister (the source image),
    %     and any additional elements are paths to other images to apply
    %     the transform to
    % baseline_nuf : char/str, optional
    %     Path to the baseline MRI to coregister to. If empty, the
    %     baseline MRI is looked up in the subject directory.
    % prefix : char/str, optional
    %     Prefix to prepend to the output filenames. Empty by default
    %     (infiles are overwritten)
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is true
    % ------------------------------------------------------------------
    arguments
        infiles
        baseline_nuf {mustBeText} = ''
        prefix {mustBeText} = ''
        overwrite logical = true
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);

    % If infiles are already at baseline, return
    mri_dir = abspath(fileparts(infiles{1}));
    subj_dir = fileparts(mri_dir);
    if isempty(baseline_nuf)
        baseline_nuf = get_baseline_mri(subj_dir);
    end
    mustBeFile(baseline_nuf);
    if strcmp(infiles{1}, baseline_nuf)
        fprintf('- This is the baseline MRI\n')
        outfiles = infiles_cp;
        return
    end

    % If outfiles already exist and overwrite is false, return
    outfiles = add_presuf(infiles, prefix);
    if all(isfile(outfiles)) && ~overwrite
        fprintf('- Will not coregister MRI to baseline, as output files already exist\n')
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        fprintf('- Coregistering MRI files to baseline MRI\n');
        for ii = 1:numel(infiles)
            if ii == 1
                fprintf('  * Source image: %s\n', basename(infiles{ii}));
            elseif ii == 2
                fprintf('  * Other images: %s\n', basename(infiles{ii}));
            else
                fprintf('                  %s\n', basename(infiles{ii}));
            end
        end
    end

    % Copy the input files if a prefix is specified
    if ~isempty(prefix)
        for ii = 1:numel(infiles)
            copyfile(infiles{ii}, outfiles{ii});
        end
    end



    % Coregister images to the baseline MRI
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(baseline_nuf);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(outfiles{1});
    if numel(outfiles) > 1
        matlabbatch{1}.spm.spatial.coreg.estimate.other = outfiles(2:end);
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run', matlabbatch);

    % Format the output filenames
    outfiles = format_outfiles(infiles_cp, prefix);
end
