function reset_origin_axis_midpoint(infiles, verbose, prefix)
    % Reset origin of each image to the midpoint along each axis
    %
    % Only the affine transform in the nifti headers are changed; data
    % arrays are unaffected.
    %
    % Parameters
    % ----------
    % infiles : str/char or cell array of str/chars of nifti filenames
    %     One or more .nii images to reorient
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

    % Reset origin of the affine transform in each image header
    for ii = 1:length(outfiles)
        % Load the image volume
        hdr = spm_vol(outfiles{ii});

        % Compute the inverse of the affine transformation matrix
        aff_inv = hdr.mat \ eye(4);

        % Update translation to the midpoint of the image dimensions
        aff_inv(1:3,4) = (hdr.dim + 1) / 2;

        % Update the affine transformation using the modified matrix
        aff = inv(aff_inv);
        spm_get_space(hdr.fname, aff);

        % Print diagnostic information
        if verbose
            fprintf('- Reset origin to axis midpoint: %s\n', basename(outfiles{ii}));
        end
    end
end
