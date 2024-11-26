function outfiles = save_pet_suvrs(pet_file, mri_dir, fid, overwrite, ref_region_file)
    % Save SUVRs in native MRI space
    %
    % Parameters
    % ----------
    % pet_file : char|string
    %     Path to the PET image that will be used for intensity
    %     normalization
    % mri_dir : char|string
    %     Path to the directory containing the MRI data
    % fid : int, optional
    %     File identifier for logging. Default is 1 for stdout
    % overwrite : logical
    %     If true, overwrite existing files
    % ref_region_file : char|string
    %     CSV file containing a mapping between PET tracer and one or
    %     more reference regions that will be used. Columns must be
    %     "tracer", "ref_region", and "masks".
    %
    % Returns
    % -------
    % outfiles : struct
    %   Struct array with paths to each output file
    % ------------------------------------------------------------------
    arguments
        pet_file {mustBeFile}
        mri_dir {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        ref_region_file {mustBeText} = ''
    end

    % Format paths
    pet_file = abspath(pet_file);
    pet_dir = fileparts(pet_file);
    pet_tag = get_scan_tag(pet_file);

    % Find the MRI directory
    if isempty(mri_dir)
        mri_dir = fullfile(pet_dir, 'mri');
    end
    mri_tag = get_scan_tag(mri_dir);

    % Get the list of reference regions and the output SUVR filenames
    [ref_regions, suvr_files] = get_suvr_files(pet_dir, ref_region_file);

    % Define the output filenames
    outfiles.rrmeans = fullfile(pet_dir, append('r', pet_tag, '_ref-region-means.csv'));
    outfiles = catstruct(outfiles, suvr_files);

    % Check if output files already exist
    if ~overwrite
        if all(structfun(@(x) isfile(x), outfiles))
            log_append(fid, '- Native MRI space SUVR files exist, will not overwrite');
            return
        end
    end

    % Initialize a table that will store outfiles.rrmeans values once we
    % have calculated them
    nrow = height(ref_regions);
    ncol = 4;
    rrmeans = table( ...
        'Size', [nrow, ncol], ...
        'VariableTypes', {'string', 'string', 'double', 'string'}, ...
         'VariableNames', {'image_file', 'mask_file', 'mean', 'voxel_count'} ...
    );

    % Iterate over rows of ref_regions. For each row, calculate the mean
    % PET value within the contributing mask(s) that define a reference
    % region and save an SUVR. Also store the mean PET value used to
    % normalize the PET image in the rrmeans table. Note that if more
    % than one mask is used to define a reference region, the reference
    % region mean is the unweighted mean of the means across the masks.
    log_append(fid, '- Saving native MRI space SUVR images');
    for ii = 1:height(ref_regions)
        % Define the output image that will be created
        ref_region = ref_regions.ref_region{ii};

        % Get path(s) to 1+ mask files that will be used to calculate
        % the reference region mean
        masks = ref_regions.masks{ii};
        maskfs = cellfun( ...
            @(x) fullfile(mri_dir, append(mri_tag, '_mask-', x, '.nii')), ...
            strtrim(split(masks, ';')), ...
            'UniformOutput', false ...
        );

        % Calculate the mean PET value within the reference region
        mask_means = roi_desc(pet_file, maskfs);

        % Convert the voxel_count column from double to string
        mask_means.voxel_count = arrayfun( ...
            @(x) sprintf('%d', x), mask_means.voxel_count, ...
            'UniformOutput', false ...
        );

        % Assign values to the rrmeans table
        rrmeans.image_file(ii) = pet_file;
        rrmeans.mask_file(ii) = strjoin(maskfs, ';\n');
        rrmeans.voxel_count(ii) = strjoin(mask_means.voxel_count, ';\n');

        % If more than one mask is used to define the reference region,
        % calculate the mean of the means (i.e. the unweighted mean
        % across the masks, regardless of their individual volumes)
        if height(mask_means) > 1
            unweighted_mean = mean(mask_means.mean);
            rrmeans.mean(ii) = unweighted_mean;
        else
            rrmeans.mean(ii) = mask_means.mean;
        end

        % Divide voxelwise PET values by the reference region mean and save the
        % output SUVR image
        ref_region = strrep(ref_region, '-', '_');
        nii_divide_img_by_val( ...
            pet_file, rrmeans.mean(ii), outfiles.(['suvr_' ref_region]), fid, overwrite ...
        );
    end

    % Save the rrmeans table to a CSV file
    writetable(rrmeans, outfiles.rrmeans);
    log_append(fid, sprintf('  * Saved %s', basename(outfiles.rrmeans)));
end
