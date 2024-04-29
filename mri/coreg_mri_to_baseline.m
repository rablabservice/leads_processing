function coreg_mri_to_baseline(mri_files, verbose)
    % Coregister MRI to baseline and ovewrite input image headers
    %
    % First looks up if mri_files already correspond to the baseline
    % MRI, in which case nothing is done.
    %
    % Parameters
    % ----------
    % mri_files : char/str or cell array of char/str
    %     Cell array of paths to the MRI files to coregister. The first
    %     element should be the MRI to coregister (the source image),
    %     and any additional elements are paths to other images to apply
    %     the transform to.
    % ------------------------------------------------------------------
    arguments
        mri_files {mustBeText}
        verbose logical = true
    end

    % Make sure all mri_files exist
    mustBeFile(mri_files);

    % If mri_files are already at baseline, return
    mri_dir = abspath(fileparts(mri_files{1}));
    subj_dir = fileparts(mri_dir);
    baseline_nuf = get_baseline_mri(subj_dir);
    mustBeFile(baseline_nuf);
    if strcmp(mri_files{1}, baseline_nuf)
        return
    end

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


    % Coregister MRI to the baseline MRI
