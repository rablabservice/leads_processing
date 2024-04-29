function outfiles = reset_origin_mri_com(infiles, verbose, prefix)
    % Recenter T1 to center-of-mass then coregister to the SPM template
    %
    % Strictly rigid-body coregistration (no reslicing is done). Only
    % the affine transform in the nifti headers are changed; data arrays
    % are unaffected.
    %
    % Parameters
    % ----------
    % infiles : str/char or cell array of str/chars of nifti filenames
    %     Images to reorient, assumed to be from same individual and
    %     session. The first image should be the T1 MRI, and the
    %     transform is estimated only for this image. Same transform is
    %     applied to any additional images.
    % verbose : logical, optional
    %     If true, print diagnostic information
    % prefix : str/char, optional
    %     Prefix to prepend to the output filenames. Empty by default
    %     (infiles are overwritten).
    % ------------------------------------------------------------------
    arguments
        infiles {mustBeText}
        verbose logical = true
        prefix {mustBeText} = ''
    end

    % Ensure infiles is a flattened cell array
    if ~iscell(infiles)
        infiles = cellstr(infiles);
    end
    infiles = infiles(:);

    % Check that all input files exist, and format them correctly
    for ii = 1:length(infiles)
        infiles{ii} = abspath(infiles{ii});
        if ~exist(infiles{ii}, 'file')
            error('File not found: %s', infiles{ii});
        end
    end

    % Print diagnostic information
    if verbose
        fprintf('- Resetting origin to center-of-mass\n')
        for ii = 1:size(infiles)
            if ii == 1
                fprintf('  - Source image: %s\n', basename(infiles{ii}));
            elseif ii == 2
                fprintf('    Other images: %s\n', basename(infiles{ii}));
            else
                fprintf('                  %s\n', basename(infiles{ii}));
            end
        end
    end

    % Copy the input files if a prefix is specified
    outfiles = cell(size(infiles));
    for ii = 1:length(infiles)
        if isempty(prefix)
            outfiles{ii} = infiles{ii};
        else
            [pth, nam, ext] = fileparts(infiles{ii});
            outfiles{ii} = fullfile(pth, [prefix nam ext]);
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
    disp(XYZ_mm);

    % Update the origin in the header of each image
    for ii = 1:length(outfiles)
        fname = outfiles{ii};
        hdr = spm_vol(fname);
        hdr.mat(1:3,4) = hdr.mat(1:3,4) - XYZ_mm(1:3);
        spm_create_vol(hdr);
    end

    % Coregister images to the SPM template
    if verbose
        fprintf('- Coregistering %s to the SPM12 OldNorm T1 template\n', basename(outfiles{1}));
    end
    coreg_to_t1(outfiles);
end


function coreg_to_t1(infiles)
    % Coregister images to the SPM template
    spm_jobman('initcfg');
    template = fullfile(spm('Dir'), 'toolbox', 'OldNorm', 'T1.nii');
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = infiles(1);
    matlabbatch{1}.spm.spatial.coreg.estimate.other = infiles(2:end);
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run', matlabbatch);
end
