function mask = create_nii_mask(infile, mask_idx, outfile, overwrite, verbose)
    % Create a binary mask from a nifti file and array of integer labels
    %
    % Mask is 1 for all elements in the input file whose values are in
    % mask_idx, and 0 otherwise.
    %
    % Parameters
    % ----------
    % infile : str|char
    %     Path to the input nifti file.
    % mask_idx : array
    %     Array of integers that represent the indices of the elements
    %     that should be included in the mask.
    % outfile : char
    %     Path to the output nifti file.
    % overwrite : bool
    %     If true, overwrite the output file if it already exists.
    % verbose : bool
    %     If true, print status messages.
    %
    % Returns
    % ------------------------------------------------------------------
    arguments
        infile {mustBeFile}
        mask_idx {mustBeNumeric}
        outfile {mustBeText} = ''
        overwrite logical = false
        verbose logical = true
    end

    % If the output file exists and overwrite is false, load the outfile
    % and return its data array
    if exist(outfile, 'file') && ~overwrite
        if verbose
            fprintf('  - File already exists, will not overwrite: %s\n', basename(outfile));
        end
        mask = spm_read_vols(spm_vol(outfile));
        return
    end

    % Make sure mask_idx is not empty
    if isempty(mask_idx)
        error('mask_idx cannot be empty');
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
    mask = ismember(dat, mask_idx);

    % Save the mask
    if save_output
        mask_img = img;
        mask_img.fname = outfile;
        mask_img.dt = [spm_type('uint8') 0]; % uint8 type, native format
        spm_write_vol(mask_img, mask);
        if verbose
            fprintf('  - Saved %s\n', basename(outfile));
        end
    end
end
