function outfiles = save_mask_pons(brainstemf, fid, overwrite, in_dir, out_dir)
    % Load the brainstem sublabels file and save the pons
    %
    % Parameters
    % ----------
    % brainstemf : char or str array
    %   Path to the brainstem-sublabels.nii file
    % fid : int, optional
    %   File identifier for logging (default is 1 for stdout)
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % in_dir : char or str array
    %   The input directory. If brainstemf is empty, this is where the
    %   function looks for the brainstem-sublabels.nii file. This
    %   parameter is disregarded if brainstemf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the brainstem-sublabels.nii file
    %
    % Files created
    % -------------
    % - <scan_tag>_mask-pons.nii
    % ------------------------------------------------------------------
    arguments
        brainstemf {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Define brainstem sublabel index for the pons
    mask_idx = 174;

    % Format inputs
    [brainstemf, out_dir] = format_mask_inputs(brainstemf);
    scan_tag = get_scan_tag(brainstemf);

    % Save the mask
    outfiles.mask_pons = fullfile(out_dir, append(scan_tag, '_mask-pons.nii'));
    nii_labels_to_mask(brainstemf, mask_idx, outfiles.mask_pons, fid, overwrite);
end
