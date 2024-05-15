function outfiles = save_mask_brainstem(aparcf, overwrite, in_dir, out_dir)
    % Load the aparc+aseg and save the brainstem mask
    %
    % Parameters
    % ----------
    % aparcf : char or str array
    %   Path to the aparc+aseg.nii file
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % in_dir : char or str array
    %   The input directory. If aparcf is empty, this is where the
    %   function looks for the aparc+aseg.nii file. This parameter is
    %   disregarded if aparcf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the aparc+aseg.nii file
    %
    % Files created
    % -------------
    % - <scan_tag>_mask-brainstem.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        overwrite logical = false
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Define aparc index for the brainstem
    mask_idx = 16;

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Save the mask
    outfiles.mask_bstem = fullfile(out_dir, append(scan_tag, '_mask-brainstem.nii'));
    nii_labels_to_mask(aparcf, mask_idx, outfiles.mask_bstem, overwrite);
end
