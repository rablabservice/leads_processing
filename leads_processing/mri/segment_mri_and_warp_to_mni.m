function segment_mri_and_warp_to_mni(mri_dir, overwrite, verbose)
    % Segment nu.nii and warp to MNI space, saving deformation files.
    %
    % Also smooths the mwc1 image by 8mm^3 FWHM to create the smwc1
    % image.
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %   The directory that contains the processed MRI data
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % <mri_dir>/c1<scan_tag>_nu.nii
    % <mri_dir>/c2<scan_tag>_nu.nii
    % <mri_dir>/c3<scan_tag>_nu.nii
    % <mri_dir>/w<scan_tag>_nu.nii
    % <mri_dir>/y_<scan_tag>_nu.nii
    % <mri_dir>/iy_<scan_tag>_nu.nii
    % <mri_dir>/smwc1<scan_tag>_nu_mni.nii
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);
    [~, nuf] = get_freesurfer_files(mri_dir, 'nii');

    % Check if segmentation/warp has already been run
    c1f = fullfile(mri_dir, append('c1', basename(nuf)));
    wnuf = fullfile(mri_dir, append('w', basename(nuf)));
    if exist(c1f, 'file') && exist(wnuf, 'file') && ~overwrite
        return
    end

    % Segment nu.nii and warp to MNI space
    if verbose
        fprintf('- Segmenting nu.nii and warping to MNI space\n');
    end
    % DO STUFF
end
