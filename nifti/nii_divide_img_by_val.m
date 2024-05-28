function outdat = nii_divide_img_by_val(infile, val, outfile, fid, overwrite)
    % Divide a nifti data array by a scalar value and save output nifti
    %
    % Parameters
    % ----------
    % infile : char|string
    %     Path to the input nifti file
    % val : numeric
    %     Scalar value to divide the data array by
    % outfile : char
    %     Path to the output nifti file
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : bool
    %     If true, overwrite the output file if it already exists
    %
    % Returns
    % outdat : numeric array
    %     Data array of the output file
    % ------------------------------------------------------------------
    arguments
        infile {mustBeFile}
        val {mustBeNumeric}
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
    if isempty(outfile)
        save_output = false;
    else
        save_output = true;
        outfile = abspath(outfile);
    end

    % Load the infile
    img = spm_vol(infile);
    dat = spm_read_vols(img);

    % Divide the data array by the scalar value
    log_append(fid, sprintf('    - Dividing %s by %f', basename(infile), val));
    outdat = dat / val;

    % Save outdat
    if save_output
        outimg = img;
        outimg.fname = outfile;
        outimg.dt = [spm_type('float32') 0];
        spm_write_vol(outimg, outdat);
        log_append(fid, sprintf('    - Saved %s', basename(outfile)));
    end
end
