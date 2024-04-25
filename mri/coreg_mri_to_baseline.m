function coreg_mri_to_baseline(mri_dir, segment_brainstem)
    % Run a single MRI scan through all the processing steps.
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        segment_brainstem logical = true
    end

    % Format inputs
    mri_dir = abspath(mri_dir);
    subj_dir = fileparts(mri_dir);

    % Get the source files
    if segment_brainstem
        [~, nuf, aparcf, brainstemf] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
        mri_files = {nuf, aparcf, brainstemf};
    else
        [~, nuf, aparcf] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
        mri_files = {nuf, aparcf};
    end

    % Get the baseline MRI and make sure it exists
    baseline_nuf = get_baseline_mri(subj_dir);
    mustBeFile(baseline_nuf);

    % Coregister MRI to the baseline MRI
