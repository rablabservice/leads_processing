function outfile = calculate_centaur(wsuvr_file,...
    outfile, ...
    fid, ...
    overwrite...
    ) 
    % Calculate CenTauR z-score values from tau PET warped images.
    %
    % First calculate SUVR in CenTauR reference region and extract
    % SUVR values from all target regions. Then apply the appropriate
    % CenTauRz conversion equation from the CenTauR manuscript.
    % For reference see Villemagne et al Alz&Dem 2023 [PMID:37424964]
    % ------------------------------------------------------------------
    arguments
        wsuvr_file
        outfile {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    code_dir = fileparts(fileparts(mfilename('fullpath')));

    % Get scan info
    wsuvr_file = abspath(cellvec(wsuvr_file));
    pet_dir = fileparts(wsuvr_file);
    pet_tag = get_scan_tag(char(pet_dir));
    [~, tracer] = parse_scan_tag(pet_tag);

    % Get the output filename and check if it exists
    if isempty(outfile)
        outfile = fullfile(pet_dir, append('wr', pet_tag, '_tau-centaur.csv'));
    end
    if isfile(outfile) && ~overwrite
        log_append(fid, sprintf('  * %s exists, will not overwrite', basename(outfile)));
        return
    end

    log_append(fid, '  * Calculating SUVRs and CenTauRz');
    
    % Get the reference region mask and value
    ref_region_maskf = fullfile(code_dir, 'centaur/roi_resliced', 'rvoi_CerebGry_tau_2mm.nii');
    ref_value = roi_desc(wsuvr_file, ref_region_maskf);
    
    % Get the target region mask and values
    frontal_maskf = fullfile(code_dir, 'centaur/roi_resliced', 'rFrontal_CenTauR.nii');
    mesial_maskf = fullfile(code_dir, 'centaur/roi_resliced', 'rMesial_CenTauR.nii');
    meta_maskf = fullfile(code_dir, 'centaur/roi_resliced', 'rMeta_CenTauR.nii');
    tp_masf = fullfile(code_dir, 'centaur/roi_resliced', 'rTP_CenTauR.nii');
    universal_maskf = fullfile(code_dir, 'centaur/roi_resliced', 'rCenTauR.nii');

    frontal_value = roi_desc(wsuvr_file, frontal_maskf);
    mesial_value = roi_desc(wsuvr_file, mesial_maskf);
    meta_value = roi_desc(wsuvr_file, meta_maskf);
    tp_value = roi_desc(wsuvr_file, tp_masf);
    universal_value = roi_desc(wsuvr_file, universal_maskf);

    % Get the intercepat and slope values
    parameters = readtable(fullfile(code_dir, 'centaur', 'centaur_constants.csv'));
    parameters = parameters(strcmpi(parameters.Tracer, tracer), :);

    % Calculate SUVR values
    frontal_suvr = frontal_value.mean/ref_value.mean;
    mesial_suvr = mesial_value.mean/ref_value.mean;
    meta_suvr = meta_value.mean/ref_value.mean;
    tp_suvr = tp_value.mean/ref_value.mean;
    universal_suvr = universal_value.mean/ref_value.mean;

    % Calculate CenTauRz values
    frontal_centaurz = parameters.Frontal_inter + (parameters.Frontal_slope * frontal_suvr);
    mesial_centaurz = parameters.MesialTemporal_inter + (parameters.MesialTemporal_slope * mesial_suvr);
    meta_centaurz = parameters.MetaTemporal_inter + (parameters.MetaTemporal_slope * meta_suvr);
    tp_centaurz = parameters.TemporoParietal_inter + (parameters.TemporoParietal_slope * tp_suvr);
    universal_centaurz = parameters.Universal_inter + (parameters.Universal_slope * universal_suvr);

    % Combine values into a table
    filenames = {wsuvr_file, wsuvr_file, wsuvr_file, wsuvr_file, wsuvr_file}; 
    masks = {universal_maskf, meta_maskf, frontal_maskf, mesial_maskf, tp_masf};     
    regions = ["Universal"; "MetaTemporal"; "Frontal"; "MesialTemporal"; "TemporoParietal"];
    suvr_values = [universal_suvr; meta_suvr; frontal_suvr; mesial_suvr; tp_suvr]; 
    centaurz_values = [universal_centaurz; meta_centaurz; frontal_centaurz; mesial_centaurz; tp_centaurz];

    % Create the table
    centaur_table = table(...
        filenames(:), masks(:), regions(:), suvr_values(:), centaurz_values(:), ...
        'VariableNames', {'Filename', 'Mask', 'Region', 'SUVR_ref-centaur', 'CenTauRz'});

    % Save the output CSV
    writetable(centaur_table, outfile);
    log_append(fid, sprintf('  * Saved %s', basename(outfile)));
end








    
