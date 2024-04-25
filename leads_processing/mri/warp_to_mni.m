function warp_to_mni(scan, mri_dir, overwrite, verbose)
    % Warp any file to MNI space based on existing y_file
    %
    %
    % Parameters
    % ----------
    % scan : char or str array
    %   The scan that needs to be warped to MNI space
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % <mri_dir>/w<scan_tag>_nu.nii
    % ------------------------------------------------------------------
    arguments
        scan {mustBeFile}
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(scan);

    % Check if warp has already been run
    [~, nuf, ~] = get_freesurfer_files(mri_dir, 'nii', 0);

    iy = fullfile(mri_dir, append('iy_', basename(nuf)));
    wf = fullfile(mri_dir, append('w', basename(scan)));
    if exist(wf, 'file') && ~overwrite
        return
    end
    
    % Check if iy_file exist
    if ~exist(iy,"file")
        error('No y-file exist! Rung segmentation first.')
    end

    if verbose
        fprintf('- Warping %s to MNI space\n', scan_tag);
    end
    
    % apply pushforward transformation using iy_file (subject -> MNI)
    % write in space of MNI152 template
    clear matlabbatch;spm('defaults','PET');
    matlabbatch{1}.spm.util.defs.comp{1}.def = cellstr(iy);
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames = cellstr(scan);
    matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = cellstr(mri_dir);
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {'/mnt/coredata/Projects/Resources/templates/icbm152.nii'};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
    spm_jobman('run',matlabbatch);
end



    