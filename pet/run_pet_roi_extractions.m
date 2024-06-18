function outfiles = run_pet_roi_extractions(suvr_files, maskfs, aparcf, roif, fid, overwrite)
    % Extract ROI means from one or more PET SUVR images
    % ------------------------------------------------------------------
    arguments
        suvr_files
        maskfs = {}
        aparcf = ''
        roif = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % Define path to the extract_rois.py program
    code_dir = fileparts(fileparts(mfilename('fullpath')));
    extract_rois_py = fullfile(code_dir, 'nifti', 'extract_rois.py');

    % Get the output file paths
    suvr_fields = fieldnames(suvr_files);
    roi_fields = cellfun( ...
        @(x) strrep(x, 'suvr', 'roi_extractions'), ...
        suvr_fields, ...
        'UniformOutput', false ...
    );
    for ii = 1:length(suvr_fields)
        suvr_field = suvr_fields{ii};
        roi_field = roi_fields{ii};
        outfiles.(roi_field) = strrep( ...
            suvr_files.(suvr_field), '.nii', '_roi-extractions.csv' ...
        );
    end

    % Check if the outputs exist
    if all(isfile(cellvec(outfiles))) && ~overwrite
        log_append(fid, '- ROI extraction files exist, will not overwrite');
        return
    end

    % Format paths
    suvr_files = abspath(cellvec(suvr_files));
    pet_dir = fileparts(suvr_files{1});
    pet_tag = get_scan_tag(pet_dir);
    [~, tracer] = parse_scan_tag(pet_tag);

    if ~isempty(maskfs)
        maskfs = abspath(cellvec(maskfs));
    end

    if ~isempty(aparcf)
        aparcf = abspath(aparcf);
        % Get the aparc ROI labels file if it wasn't provided
        if isempty(roif)
            roif = fullfile(code_dir, 'config', sprintf('fsroi_list_%s.csv', tracer));
        end
    end

    % Extract ROI means from each PET image
    log_append(fid, '- Extracting ROI means from PET SUVRs in native MRI space:');
    for ii = 1:length(suvr_files)
        % Construct the ROI extraction command
        suvr_file = suvr_files{ii};
        roi_field = roi_fields{ii};
        outfile = outfiles.(roi_field);
        cmd = sprintf('%s --images %s', extract_rois_py, suvr_file);
        if ~isempty(maskfs)
            cmd = sprintf('%s --masks %s', cmd, strjoin(maskfs));
        end
        if ~isempty(aparcf)
            cmd = sprintf('%s --aparcs %s --roi_file %s', cmd, aparcf, roif);
        end
        cmd = sprintf('%s --shape long --outputf %s --quiet', cmd, outfile);
        fprintf(cmd);

        % Run the command
        system(cmd);
    end
end
