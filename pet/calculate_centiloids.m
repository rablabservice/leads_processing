function outfile = calculate_centiloids(suvr_files, cortical_summary_maskf, outfile, fid, overwrite)
    % Calculate Centiloid values from amyloid PET SUVR images.
    %
    % First calculate SUVR means in native space within the ADNI
    % cortical summary regions (comprised of a set of FreeSurfer ROI
    % labels from the Desikan-Killiany atlas). Then apply the
    % appropriate Centiloid conversion equation from Royce et al. 2021:
    % doi.org/10.1186/s13195-021-00836-1
    % ------------------------------------------------------------------
    arguments
        suvr_files
        cortical_summary_maskf {mustBeText} = ''
        outfile {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % Get scan info
    suvr_files = abspath(cellvec(suvr_files));
    pet_dir = fileparts(suvr_files{1});
    pet_tag = get_scan_tag(pet_dir);
    [~, tracer] = parse_scan_tag(pet_tag);

    % Get the output filename and check if it exists
    if isempty(outfile)
        outfile = fullfile(pet_dir, append('r', pet_tag, '_amyloid-cortical-summary.csv'));
    end
    if isfile(outfile) && ~overwrite
        log_append(fid, sprintf('  * %s exists, will not overwrite', basename(outfile)));
        return
    end

    % Get the cortical summary mask if not provided
    if isempty(cortical_summary_maskf)
        mri_dir = fullfile(pet_dir, 'mri');
        mri_tag = get_scan_tag(mri_dir);
        cortical_summary_maskf = fullfile( ...
            pet_dir, 'mri', append(mri_tag, '_mask-', 'amyloid-cortical-summary.nii') ...
        );
    end

    % Extract SUVR means from the cortical summary mask
    log_append(fid, '  * Calculating cortical summary SUVRs and Centiloids');
    cl_table = roi_desc(suvr_files, cortical_summary_maskf);
    cl_table.Properties.VariableNames{'mean'} = 'mean_suvr';
    cl_table.voxel_count = [];
    cl_table.centiloids = NaN(height(cl_table), 1);

    % Apply the Centiloid conversion equations from Royce et al. 2021
    if strcmp(tracer, 'FBB')
        idx = contains(cl_table.image_file, 'suvr-wcbl');
        cl_table(idx, :).centiloids = (157.15 * cl_table(idx, :).mean_suvr) - 151.87;
        idx = contains(cl_table.image_file, 'suvr-compwm');
        cl_table(idx, :).centiloids = (244.20 * cl_table(idx, :).mean_suvr) - 170.80;
    elseif strcmp(tracer, 'FBP')
        idx = contains(cl_table.image_file, 'suvr-wcbl');
        cl_table(idx, :).centiloids = (188.22 * cl_table(idx, :).mean_suvr) - 189.16;
        idx = contains(cl_table.image_file, 'suvr-compwm');
        cl_table(idx, :).centiloids = (300.66 * cl_table(idx, :).mean_suvr) - 208.84;
    else
        warning('No Centiloid conversion equation has been programmed for %s', tracer);
    end

    % Save the output CSV
    writetable(cl_table, outfile);
    log_append(fid, sprintf('  * Saved %s', outfile));
end
