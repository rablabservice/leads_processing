function mrifs = get_processed_mri_files(mri_dir)
    % Check which processed MRI files exist for a given MRI directory
    %
    % Parameters
    % ----------
    % mri_dir : char or str
    %     The directory that contains the processed MRI data
    %
    % Returns
    % -------
    % mrifs_exist : struct
    %     A struct containing the paths to the processed MRI files that
    %     exist
    % mrifs_missing : struct
    %     A struct containing the paths to the processed MRI files that
    %     do not exist
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText} = ''
    end

    % Get scan info
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % Define the processed MRI files
    mrifs = struct( ...
        'nu', fullfile( ...
            mri_dir, append(scan_tag, '_nu.nii') ...
        ), ...
        'aparc', fullfile( ...
            mri_dir, append(scan_tag, '_aparc+aseg.nii') ...
        ), ...
        'bstem', fullfile( ...
            mri_dir, append(scan_tag, '_brainstem_sublabels.nii') ...
        ), ...
        'suit_atlas', fullfile( ...
            mri_dir, append(scan_tag, '_cbl-suit.nii') ...
        ), ...
        'mask_wcbl', fullfile( ...
            mri_dir, append(scan_tag, '_mask-wcbl.nii') ...
        ), ...
        'mask_bstem', fullfile( ...
            mri_dir, append(scan_tag, '_mask-brainstem.nii') ...
        ), ...
        'mask_cortical_summary', fullfile( ...
            mri_dir, append(scan_tag, '_mask-amyloid-cortical-summary.nii') ...
        ), ...
        'mask_subcortwm', fullfile( ...
            mri_dir, append(scan_tag, '_mask-subcortwm.nii') ...
        ), ...
        'ssubcortwm', fullfile( ...
            mri_dir, append(scan_tag, '_mask-subcortwm.nii') ...
        ), ...
        'mask_eroded_subcortwm', fullfile( ...
            mri_dir, append(scan_tag, '_mask-eroded-subcortwm.nii') ...
        ), ...
        'mask_cblgm', fullfile( ...
            mri_dir, append(scan_tag, '_mask-cblgm.nii') ...
        ), ...
        'mask_infcblgm', fullfile( ...
            mri_dir, append(scan_tag, '_mask-infcblgm.nii') ...
        ), ...
        'pons', fullfile( ...
            mri_dir, append(scan_tag, '_mask-pons.nii') ...
        ), ...
        'c1', fullfile( ...
            mri_dir, append('c1', scan_tag, '_nu.nii') ...
        ), ...
        'c2', fullfile( ...
            mri_dir, append('c2', scan_tag, '_nu.nii') ...
        ), ...
        'c3', fullfile( ...
            mri_dir, append('c3', scan_tag, '_nu.nii') ...
        ), ...
        'c4', fullfile( ...
            mri_dir, append('c4', scan_tag, '_nu.nii') ...
        ), ...
        'c5', fullfile( ...
            mri_dir, append('c5', scan_tag, '_nu.nii') ...
        ), ...
        'iy', fullfile( ...
            mri_dir, append('iy_', scan_tag, '_nu.nii') ...
        ), ...
        'smwc1', fullfile( ...
            mri_dir, append('smwc1', scan_tag, '_nu.nii') ...
        ), ...
        'wnu', fullfile( ...
            mri_dir, append('w', scan_tag, '_nu.nii') ...
        ), ...
        'y', fullfile( ...
            mri_dir, append('y_', scan_tag, '_nu.nii') ...
        ), ...
        'atf', fullfile( ...
            mri_dir, append('atf_', scan_tag, '_nu.mat') ...
        ), ...
        'anu', fullfile( ...
            mri_dir, append('a', scan_tag, '_nu.nii') ...
        ) ...
    );
end