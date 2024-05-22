function outfiles = save_mask_wcbl(aparcf, fid, overwrite, in_dir, out_dir)
    % Load the aparc+aseg and save the whole cerebellum mask
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
    % - <scan_tag>_mask-wcbl.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Define aparc indices for the whole cerebellum
    mask_idx = [7; 8; 46; 47];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Save the mask
    outfiles.mask_wcbl = fullfile(out_dir, append(scan_tag, '_mask-wcbl.nii'));
    nii_labels_to_mask(aparcf, mask_idx, outfiles.mask_wcbl, fid, overwrite);
end
