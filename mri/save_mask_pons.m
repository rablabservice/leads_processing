function outfile = save_mask_pons(brainstemf, in_dir, out_dir, overwrite, verbose)
    % Load the brainstem sublabels file and save the pons
    %
    % Usage
    % -----
    % >> save_mask_pons(brainstemf, in_dir, out_dir, overwrite, verbose)
    %
    % Parameters
    % ----------
    % brainstemf : char or str array
    %   Path to the brainstem-sublabels.nii file
    % in_dir : char or str array
    %   The input directory. If brainstemf is empty, this is where the
    %   function looks for the brainstem-sublabels.nii file. This
    %   parameter is disregarded if brainstemf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the brainstem-sublabels.nii file
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % - <out_dir>/<scan_tag>_mask-pons.nii
    % ------------------------------------------------------------------
    arguments
        brainstemf {mustBeText} = ''
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
        overwrite logical = false
        verbose logical = true
    end

    % Define brainstem sublabel index for the pons
    mask_idx = 174;

    % Format inputs
    [brainstemf, out_dir] = format_mask_inputs(brainstemf);
    scan_tag = get_scan_tag(brainstemf);

    % Save the mask
    outfile = fullfile(out_dir, append(scan_tag, '_mask-pons.nii'));
    nii_labels_to_mask(brainstemf, mask_idx, outfile, overwrite, verbose);
end
