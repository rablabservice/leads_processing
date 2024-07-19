function outdat = nii_apply_mask(datfile, maskfile, invert_mask, outfile, fid, overwrite)
    % Apply mask to a nifti image. By default, datfile values outside
    % the mask are set to 0.
    %
    % Parameters
    % ----------
    % datfile : char|string|numeric array
    %     Path to the input nifti file with data that we want to mask,
    %     or the data array itself. If datfile is a numeric array, the
    %     output image will not be saved.
    % maskfile : char|string|numeric array
    %     Path to the input nifti file with the mask that we will apply
    %     to datfile, or the mask array itself.
    % invert_mask : bool, default = false
    %     If false, datfile values outside the mask are set to 0.
    %     If true, datfile values within the mask are set to 0.
    % outfile : char
    %     Path to the output nifti file. If empty, the output will not
    %     be saved unless overwrite is true, in which case datfile will
    %     be overwritten.
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : bool
    %     If true, overwrite the output file if it already exists.
    %
    % Returns
    % -------
    % outdat : numeric array
    %     Data array of the output file
    % ------------------------------------------------------------------
    arguments
        datfile
        maskfile
        invert_mask logical = false
        outfile {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % If the output file exists and overwrite is false, load the outfile
    % and return its data array
    if isfile(outfile) && ~overwrite
        log_append(fid, sprintf('    - %s exists, will not overwrite', basename(outfile)));
        outdat = spm_read_vols(spm_vol(outfile));
        return
    end

    % Determine if we should save the output
    if isempty(outfile) || isnumeric(datfile)
        save_output = false;
    else
        save_output = true;
        outfile = abspath(outfile);
    end

    % Load the data
    if isnumeric(datfile) || islogical(datfile)
        dat = datfile;
    else
        img = spm_vol(datfile);
        dat = spm_read_vols(img);
    end

    % Load the mask
    if isnumeric(maskfile) || islogical(maskfile)
        mask = maskfile;
    else
        mask = spm_read_vols(spm_vol(maskfile));
    end

    % Invert the mask, if relevant
    if invert_mask
        mask = ~mask;
    end

    % Apply the mask
    outdat = dat;
    outdat(~mask) = 0;

    % Save outdat
    if save_output
        outimg = img;
        outimg.fname = outfile;
        spm_write_vol(outimg, outdat);
        log_append(fid, sprintf('    - Saved %s', basename(outfile)));
    end
end
