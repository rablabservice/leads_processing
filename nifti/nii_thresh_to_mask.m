function mask = nii_thresh_to_mask(infile, lower, upper, outfile, fid, overwrite)
    % Create a binary mask of infile values > lower and < upper
    %
    % Mask is 1 for all elements in the input file whose values are in
    % labels, and 0 otherwise.
    %
    % Parameters
    % ----------
    % infile : char|string
    %     Path to the input nifti file
    % lower : float
    %     Threshold at or below which voxels are excluded from the mask
    % upper : float
    %     Threshold at or above which voxels are excluded from the mask
    % outfile : char
    %     Path to the output nifti file
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : bool
    %     If true, overwrite the output file if it already exists
    %
    % Returns
    % mask : logical array
    %     Mask array in the shape of the input file
    % ------------------------------------------------------------------
    arguments
        infile {mustBeFile}
        lower {mustBeNumeric} = -Inf
        upper {mustBeNumeric} = Inf
        outfile {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % If the output file exists and overwrite is false, load the outfile
    % and return its data array
    if exist(outfile, 'file') && ~overwrite
        log_append(fid, sprintf('  * %s exists, will not overwrite', basename(outfile)));
        mask = spm_read_vols(spm_vol(outfile));
        return
    end

    % Make sure at least one threshold is specified
    if ~isfinite(lower) && ~isfinite(upper)
        error('At least one of lower or upper must be finite');
    end

    % Determine if we should save the output
    if isempty(outfile)
        save_output = false;
    else
        save_output = true;
        outfile = abspath(outfile);
    end

    % Load the infile
    img = spm_vol(infile);
    dat = spm_read_vols(img);

    % Create the mask
    mask = (dat > lower) & (dat < upper);

    % Save the mask
    if save_output
        mask_img = img;
        mask_img.fname = outfile;
        mask_img.dt = [spm_type('uint8') 0];
        spm_write_vol(mask_img, mask);
        log_append(fid, sprintf('  * Saved %s', basename(outfile)));
    end
end
