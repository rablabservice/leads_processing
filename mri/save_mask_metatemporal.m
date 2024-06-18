function outfiles = save_mask_metatemporal(aparcf, fid, overwrite, in_dir, out_dir)
    % Load the aparc+aseg and save the metatemporal ROI mask
    %
    % Parameters
    % ----------
    % aparcf : char or str array
    %   Path to the aparc+aseg.nii file
    % fid : int, optional
    %   File identifier for logging (default is 1 for stdout)
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
    % - <scan_tag>_mask-metatemporal.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Define aparc indices for the metatemporal ROI (comprises left and
    % right AMY, ERC, PHG, FSF, ITG, and MTG)
    mask_idx = [18; 54; 1006; 1007; 1009; 1015; 1016; 2006; 2007; 2009; 2015; 2016];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Save the mask
    outfiles.mask_metatemporal = fullfile(out_dir, append(scan_tag, '_mask-metatemporal.nii'));
    nii_labels_to_mask(aparcf, mask_idx, outfiles.mask_metatemporal, fid, overwrite);
end
