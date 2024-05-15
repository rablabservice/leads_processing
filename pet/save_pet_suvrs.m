function outfiles = save_pet_suvrs(petf, overwrite, mri_dir, ref_regionf, metadata_dir)
    % Save SUVRs in native MRI space
    %
    % Parameters
    % ----------
    % petf : char|string
    %     Path to the PET image that will be used for intensity
    %     normalization
    % overwrite : logical
    %     If true, overwrite existing files
    % mri_dir : char|string
    %     Path to the directory containing the MRI data
    % ref_regionf : char|string
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
        petf {mustBeFile}
        overwrite logical = false
        mri_dir {mustBeText} = ''
        ref_regionf {mustBeText} = ''
        metadata_dir {mustBeText} = '/mnt/coredata/processing/leads/metadata'
    end

    % Format paths
    petf = abspath(petf);
    pet_tag = get_scan_tag(petf);
    [subj, tracer, pet_date] = parse_scan_tag(pet_tag);

    % Find the MRI directory
    if isempty(mri_dir)
        mri_dir = fullfile(dirname(petf), 'mri');
        mri_tag = get_scan_tag(mri_dir);
    end

    % Get a list of reference regions to use
    if isempty(ref_regionf)
        ref_regionf = fullfile(metadata_dir, 'ssheets', 'ref_regions.csv');
        ref_regions = readtable(ref_regionf);
    end

    % Load the PET image
    pet_img = spm_vol(petf);
    pet_dat = spm_read_vols(pet_img);

    % Find the reference region(s) that we'll use

    outfiles.suvr = fullfile(mri_dir, append('suvr_', pet_tag, '_nu.nii'));

    % Check if output files already exist
    if all(structfun(@(x) exist(x, 'file'), outfiles)) && ~overwrite
        fprintf('- SUVRs already saved, will not rerun\n')
        return
    else
        fprintf('- Saving SUVRs for %s\n', basename(petf));
    end

    % Load the PET image
    pet = spm_vol(petf);
    pet_img = spm_read_vols(pet);

    % Load the mask
    mask = spm_vol(fullfile(mri_dir, append('mask_', mri_tag, '_nu.nii')));
    mask_img = spm_read_vols(mask);

    % Calculate the SUVR
    suvr_img = pet_img ./ mean(pet_img(mask_img > 0));

    % Save the SUVR
    suvr = pet;
    suvr.fname = outfiles.suvr;
    suvr.descrip = 'SUVR image';
    spm_write_vol(suvr, suvr_img);

    % Save the SUVR mask
    suvr_mask = mask;
    suvr_mask.fname = outfiles.suvr_mask;
    suvr_mask.descrip = 'SUVR mask';
