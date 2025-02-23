function pet_files = get_processed_pet_files(pet_dir)
    % Return a struct of all processed PET files for a selected scan
    % that should exist after processing is complete
    %
    % Parameters
    % ----------
    % pet_dir : char or str
    %     The directory that contains the processed PET data
    %
    % Returns
    % -------
    % pet_files : struct
    %     A struct containing the paths to the processed PET files
    % ------------------------------------------------------------------
    arguments
        pet_dir {mustBeText} = ''
    end

    % Hard-code a list amyloid PET tracers
    amyloid_tracers = {'FBB', 'FBP', 'FLUTE', 'NAV', 'PIB'};
    tau_tracers = {'FTP', 'PI2620', 'MK6240'};

    % Get scan info
    pet_dir = abspath(pet_dir);
    pet_tag = get_scan_tag(pet_dir);
    [~, tracer] = parse_scan_tag(pet_tag);
    tracer_is_amyloid = ismember(tracer, amyloid_tracers);
    tracer_is_tau = ismember(tracer, tau_tracers);

    % Get names of all the processed PET files
    pet_files = struct( ...
        'pet', fullfile( ...
            pet_dir, append(pet_tag, '.nii') ...
        ), ...
        'rpet', fullfile( ...
            pet_dir, append('r', pet_tag, '.nii') ...
        ), ...
        'rrmeans', fullfile( ...
            pet_dir, append('r', pet_tag, '_ref-region-means.csv') ...
        ), ...
        'qc_image', fullfile( ...
            pet_dir, append(pet_tag, '_qc.png') ...
        ) ...
    );
    [~, suvr_files] = get_suvr_files(pet_dir);
    pet_files = catstruct(pet_files, suvr_files);

    suvr_fields = fieldnames(suvr_files);
    roi_fields = cellfun( ...
        @(x) strrep(x, 'suvr', 'roi_extractions'), ...
        suvr_fields, ...
        'UniformOutput', false ...
    );
    for ii = 1:length(suvr_fields)
        suvr_field = suvr_fields{ii};
        roi_field = roi_fields{ii};
        pet_files.(roi_field) = strrep( ...
            suvr_files.(suvr_field), '.nii', '_roi-extractions.csv' ...
        );
    end

    if tracer_is_amyloid
        pet_files.cortical_summary = fullfile(pet_dir, append('r', pet_tag, '_amyloid-cortical-summary.csv'));
    elseif tracer_is_tau
        pet_files.centaur = fullfile(pet_dir, append('wr', pet_tag, '_tau-centaur.csv'));
    end

    pet_files = catstruct(pet_files, format_outfiles(suvr_files, 'w'));
    pet_files = catstruct(pet_files, format_outfiles(suvr_files, 'a'));
end
