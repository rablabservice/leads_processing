function reset_origin_axis_midpoint(infiles, prefix, overwrite)
    % Reset origin of each image to the midpoint along each axis
    %
    % Only the affine transform in the nifti headers are changed; data
    % arrays are unaffected.
    %
    % Parameters
    % ----------
    % infiles : str/char or cell array of str/chars of nifti filenames
    %     One or more .nii images to reorient
    % prefix : str/char, optional
    %     Prefix to prepend to the output filenames. Empty by default
    %     (infiles are overwritten)
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is true
    % ------------------------------------------------------------------
    arguments
        infiles
        prefix {mustBeText} = ''
        overwrite logical = true
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);

    % If outfiles already exist and overwrite is false, return
    outfiles = add_presuf(infiles, prefix);
    if all(isfile(outfiles)) && ~overwrite
        fprintf('- Will not reset origin to axis midpoint, as output files already exist\n')
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        fprintf('- Resetting origin to axis midpoint\n')
        for ii = 1:numel(infiles)
            if isempty(prefix)
                fprintf('  * %s\n', basename(infiles{ii}));
            else
                fprintf('  * %s -> %s\n', basename(infiles{ii}), basename(outfiles{ii}));
            end
        end
    end

    % Copy the input files if a prefix is specified
    if ~isempty(prefix)
        for ii = 1:numel(infiles)
            copyfile(infiles{ii}, outfiles{ii});
        end
    end

    % Reset origin of the affine transform in each image header
    for ii = 1:numel(outfiles)
        % Load the image volume
        hdr = spm_vol(outfiles{ii});

        % Compute the inverse of the affine transformation matrix
        aff_inv = hdr.mat \ eye(4);

        % Update translation to the midpoint of the image dimensions
        aff_inv(1:3,4) = (hdr.dim + 1) / 2;

        % Update the affine transformation using the modified matrix
        aff = inv(aff_inv);
        spm_get_space(hdr.fname, aff);
    end

    % Format the output filenames
    outfiles = format_outfiles(infiles_cp, prefix);
end
