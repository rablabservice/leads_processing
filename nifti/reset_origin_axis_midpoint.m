function outfiles = reset_origin_axis_midpoint(infiles, fid, overwrite, prefix)
    % Reset origin of each image to the midpoint along each axis
    %
    % Only the affine transform in the nifti headers are changed; data
    % arrays are unaffected.
    %
    % Parameters
    % ----------
    % infiles : str/char or cell array of str/chars of nifti filenames
    %     One or more .nii images to reorient
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is true
    % prefix : str/char, optional
    %     Prefix to prepend to the output filenames. Empty by default
    %     (infiles are overwritten)
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
        log_append(fid, '- Will not reset origin to axis midpoint, as output files exist');
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        log_append(fid, '- Resetting origin to axis midpoint');
        for ii = 1:numel(infiles)
            if isempty(prefix)
                log_append(fid, sprintf('  * %s', basename(infiles{ii})));
            else
                msg = sprintf( ...
                    '  * %s ->\n              %s', ...
                    basename(infiles{ii}), ...
                    basename(outfiles{ii}) ...
                );
                log_append(fid, msg);
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
