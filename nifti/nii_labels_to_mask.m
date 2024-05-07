function mask = nii_labels_to_mask(infile, labels, outfile, overwrite)
    % Create a binary mask from a nifti file and array of integer labels
    %
    % Mask is 1 for all elements in the input file whose values are in
    % labels, and 0 otherwise.
    %
    % Parameters
    % ----------
    % infile : char|string
    %     Path to the input nifti file.
    % labels : array
    %     Array of integers that represent the indices of the elements
    %     that should be included in the mask.
    % outfile : char
    %     Path to the output nifti file.
    % overwrite : bool
    %     If true, overwrite the output file if it already exists.
    %
    % Returns
    % mask : logical array
    %     Mask array in the shape of the input file
    % ------------------------------------------------------------------
    arguments
        infile {mustBeFile}
        labels {mustBeNumeric}
        outfile {mustBeText} = ''
        overwrite logical = false
    end

    % If the output file exists and overwrite is false, load the outfile
    % and return its data array
    if exist(outfile, 'file') && ~overwrite
        fprintf('  * %s already exists, will not overwrite\n', basename(outfile));
        mask = spm_read_vols(spm_vol(outfile));
        return
    end

    % Make sure labels is not empty
    if isempty(labels)
        error('labels cannot be empty');
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
    mask = ismember(dat, labels);

    % Save the mask
    if save_output
        mask_img = img;
        mask_img.fname = outfile;
        mask_img.dt = [spm_type('uint8') 0];
        spm_write_vol(mask_img, mask);
        fprintf('  * Saved %s\n', basename(outfile));
    end
end
