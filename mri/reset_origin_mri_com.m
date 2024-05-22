function outfiles = reset_origin_mri_com(infiles, fid, overwrite, prefix)
    % Recenter T1 to center-of-mass then coregister to the SPM template
    %
    % Strictly rigid-body coregistration (no reslicing is done). Only
    % the affine transform in the nifti headers are changed; data arrays
    % are unaffected
    %
    % Parameters
    % ----------
    % infiles : char|string|cellstr|struct of filenames
    %     Images to reorient, assumed to be from same individual and
    %     session. The first image should be the T1 MRI. The transform
    %     is estimated for the first image but applied to all infiles
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is true
    % prefix : str/char, optional
    %     Prefix to prepend to the output filenames. Empty by default
    %     (infiles are overwritten)
    %
    % Returns
    % -------
    % outfiles : char|string|cellstr|struct of filenames
    %     Paths to the reoriented images, in the same order and object
    %     class as the input images
    % ------------------------------------------------------------------
    arguments
        infiles
        fid {mustBeNumeric} = 1
        overwrite logical = true
        prefix {mustBeText} = ''
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);

    % If outfiles already exist and overwrite is false, return
    outfiles = add_presuf(infiles, prefix);
    if all(isfile(outfiles)) && ~overwrite
        log_append(fid, '- Will not reset MRI origin to center-of-mass, as output files exist');
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        log_append(fid, '- Resetting MRI origin to center-of-mass');
        for ii = 1:numel(infiles)
            if ii == 1
                log_append(fid, sprintf('  * Source image: %s', basename(infiles{ii})));
            elseif ii == 2
                log_append(fid, sprintf('  * Other images: %s', basename(infiles{ii})));
            else
                log_append(fid, sprintf('                  %s', basename(infiles{ii})));
            end
        end
    end

    % Copy the input files if a prefix is specified
    if ~isempty(prefix)
        for ii = 1:numel(infiles)
            copyfile(infiles{ii}, outfiles{ii});
        end
    end

    % Load the first image and find the center of mass
    hdr = spm_vol(outfiles{1});
    img = spm_read_vols(hdr);
    img = img - min(img(:));
    img(isnan(img)) = 0;
    sumTotal = sum(img(:));
    coivox = ones(4,1);
    coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal;
    coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal;
    coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal;
    XYZ_mm = hdr.mat * coivox;

    % Update the origin in the header of each image
    for ii = 1:length(outfiles)
        fname = outfiles{ii};
        hdr = spm_vol(fname);
        hdr.mat(1:3,4) = hdr.mat(1:3,4) - XYZ_mm(1:3);
        spm_create_vol(hdr);
    end

    % Find the template that we'll coregister to
    templatef = fullfile(spm('Dir'), 'toolbox', 'OldNorm', 'T1.nii');
    mustBeFile(templatef);

    % Coregister images to the SPM template
    log_append(fid, '- Coregistering MRI to the OldNorm T1 template');
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(templatef);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = outfiles(1);
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
